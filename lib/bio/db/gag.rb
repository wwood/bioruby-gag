

class Bio::DB::PileupIterator
  # Find places in this pileup that correspond to GAG errors
  # * Currently XYX in the consensus (i.e. the first and third bases are the same, and different to the middle one)
  # * There is at least 3 reads that have an insertion of base Y next to Y, and are all in the one direction. Can change this with options[:min_disagreeing_absolute]
  # * The 3 or more reads form at least a proportion of 0.1 (i.e. 10%) of all the reads at that position.  Can change this with options[:min_disagreeing_proportion]
  #
  # Returns an array of Bio::Gag objects
  #
  # When a block is given, each gag is yielded
  def gags(options={})
    min_disagreeing_proportion = options[:min_disagreeing_proportion]
    min_disagreeing_proportion ||= 0.1
    min_disagreeing_absolute = options[:min_disagreeing_absolute]
    min_disagreeing_absolute ||= 3
    
    log = Bio::Log::LoggerPlus['bio-gag']
    
    piles = []
    gags = []
    
    each do |pile|
      if piles.length < 2
        log.debug "Piles cache for this reference sequence less than length 2"
        piles = [piles, pile].flatten
        next
      elsif piles.length < 3
        log.debug "Piles cache for this reference sequence becoming full"
        piles = [piles, pile].flatten
      elsif piles[1].ref_name != pile.ref_name
        log.debug "Piles cache removed - moving to new contig"
        piles = [pile]
        next
      else
        log.debug "Piles cache regular push through"
        piles = [piles[1], piles[2], pile].flatten
      end
      log.debug "Current piles now at #{piles[0].ref_name}, #{piles.collect{|pile| "#{pile.pos}/#{pile.ref_base}"}.join(', ')}"
      
      # if not at the start/end of the contig
      first = piles[0]
      second = piles[1]
      third = piles[2] 
      
      # First and third nucleic acids must be the same
      if first.ref_base.upcase != third.ref_base.upcase
        log.debug "First and third bases are not equivalent (#{first.ref_base.upcase} and #{third.ref_base.upcase}), so not gagging"
        next
      end
      
      # can't be all one homopolymer
      if first.ref_base.upcase == second.ref_base.upcase
        log.debug "First and second bases are equivalent, so not gagging"
        next
      end
      
      # all reads that have a single insertion after the first or second position, but not both  
      inserting_reads = [first.reads, second.reads].flatten.uniq.select do |read|
        !(read.insertions[first.pos] and read.insertions[second.pos]) and
        (read.insertions[first.pos] or read.insertions[second.pos])
      end
      log.debug "Inserting reads after filtering: #{inserting_reads.inspect}"
      
      # ignore regions that aren't ever going to make it past the next filter
      if inserting_reads.length < min_disagreeing_absolute or inserting_reads.length.to_f/first.coverage < min_disagreeing_proportion
        log.debug "Insufficient disagreement at step 1, so not calling a gag"
        next
      end

      # what is the maximal base that is inserted and maximal number of directions
      direction_counts = {'+' => 0, '-' => 0}
      base_counts = {}
      inserting_reads.each do |read|
        insert = read.insertions[first.pos]
        insert ||= read.insertions[second.pos]
        insert.upcase!
        direction_counts[read.direction] += 1
        base_counts[insert] ||= 0
        base_counts[insert] += 1
      end
      log.debug "Direction counts of insertions: #{direction_counts.inspect}"
      log.debug "Base counts of insertions: #{base_counts.inspect}"
      max_direction = direction_counts['+']>direction_counts['-'] ? '+' : '-'
      max_base = base_counts.max do |a,b|
        a[1] <=> b[1]
      end[0]
      log.debug "Picking max direction #{max_direction} and max base #{max_base}"
      
      # Only accept positions that are inserting a single base
      if max_base.length > 1
        log.debug "Maximal insertion is too long, so not calling a gag"
        next
      end
      
      counted_inserts = inserting_reads.select do |read|
        insert = read.insertions[first.pos]
        insert ||= read.insertions[second.pos]
        insert.upcase!
        read.direction == max_direction and insert == max_base
      end
      log.debug "Reads counting after final filtering: #{counted_inserts.inspect}"
      
      coverage = (first.coverage+second.coverage+third.coverage).to_f / 3.0
      coverage_percent = counted_inserts.length.to_f / coverage
      log.debug "Final abundance calculations: max base #{max_base} (comparison base #{second.ref_base.upcase}) occurs #{counted_inserts.length} times compared to coverage #{coverage} (#{coverage_percent*10}%)"
      if max_base != second.ref_base.upcase or # first and second bases must be the same 
        counted_inserts.length < min_disagreeing_absolute or # require 3 bases in that maximal direction
        coverage_percent < min_disagreeing_proportion # at least 10% of reads with disagree with the consensus and agree with the gag
        log.debug "Failed final abundance cutoffs, so not calling a gag"
        next          
      end 
      
      # alright, gamut navigated. We have a match, record it
      gag = Bio::Gag.new(second.pos, piles, first.ref_name)
      gags.push gag
      yield gag if block_given?
    end
    
    return gags
  end
  
  # Given a hash containing sequence identifier => sequences, where both key and value are plain old Ruby strings, return the hash with any GAG errors in the sequences fixed
  def fix_gags(hash_of_sequence_ids_to_sequence_strings)
    log = Bio::Log::LoggerPlus['bio-gag']
    
    # Get the gags
    sequence_id_to_gags = {} #Hash of sequence id to array of gag objects
    gags do |gag|
      sequence_id_to_gags[gag.ref_name] ||= []
      sequence_id_to_gags[gag.ref_name].push gag 
    end
    log.info "Found #{sequence_id_to_gags.values.flatten.length} gag errors"
    
    # Make sure all gag errors in the pileup map to a sequence input fasta file by keeping tally
    accounted_for_seq_ids = []
    fixed_sequences = {} #Hash of sequence ids to sequences without gag errors
    hash_of_sequence_ids_to_sequence_strings.each do |seq_id, seq|
      log.debug "Now attempting to fix sequence #{seq_id}, sequence #{seq}"
      toilet = sequence_id_to_gags[seq_id]
      if toilet.nil?
        # No gag errors found in this sequence (or pessimistically the sequence wasn't in the pileup -leaving that issue to the user though)
        fixed_sequences[seq_id] = seq
      else
        # Gag error found at least once somewhere in this sequence
        # Record that this was touched in the pileup
        accounted_for_seq_ids.push seq_id
        
        # Output the fixed-up sequence
        last_gag = 0
        fixed = ''
        toilet.sort{|a,b| a.position<=>b.position}.each do |gag|
          fixed = fixed+seq[last_gag..(gag.position-1)]
          fixed = fixed+seq[(gag.position-1)..(gag.position-1)]
          last_gag = gag.position
          log.debug "After fixing gag at position #{gag.position}, fixed sequence is now #{fixed}"
        end
        fixed = fixed+seq[last_gag..(seq.length-1)]
        fixed_sequences[seq_id] = fixed
      end
    end
    
    unless accounted_for_seq_ids.length == sequence_id_to_gags.length
      log.warn "Unexpectedly found GAG errors in sequences that weren't in the sequence that are to be fixed: Found gags in #{gags.length}, but only fixed #{accounted_for_seq_ids.length}"
    end
    return fixed_sequences
  end
end

class Bio::Gag
  # The name of the reference sequence where the error was called
  attr_accessor :ref_name
  
  # Position in the reference sequence where the error was called
  attr_accessor :position
  
  # Bio::DB::Pileup objects around the GAG error
  attr_accessor :gagging_pileups
  
  def initialize(position, gagging_pileups, ref_name)
    @position = position
    @gagging_pileups = gagging_pileups
    @ref_name = ref_name
  end
  
  def inserted_base
    @gagging_pileups[1].ref_base
  end
end


