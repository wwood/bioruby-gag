

class Bio::DB::PileupIterator
  DEFAULT_GAG_ERROR_CONTEXTS = %w(GAG CTC  AGC GCT  GCG CGC  GCA TGC)
  ALL_POSSIBLE_GAG_ERROR_CONTEXTS = %w(
      GAG GTG GCG
      CTC CAC CGC
      AGA ACA ATA
      TAT TGT TCT
      )
  
  # Find places in this pileup that correspond to GAG errors
  # * Only certain sequences are considered to be possible errors. Can change this with options[:acceptable_gag_errors]
  # ** GAAG/CTTC (namesake of GAG errors. I.e. GAG is looked for, to see if it is probably GAAG instead)
  # ** AGGC/GCCT
  # ** GCCG/CGGC
  # ** GCCA/TGGC
  # * There is at least 3 reads that have an insertion of base Y next to Y, and are all in the one direction. Can change this with options[:min_disagreeing_absolute]
  # * The 3 or more reads form at least a proportion of 0.1 (i.e. 10%) of all the reads at that position.  Can change this with options[:min_disagreeing_proportion]
  #
  # To not restrict options[:acceptable_gag_errors] to any sequences, use  Bio::DB::PileupIterator::ALL_POSSIBLE_GAG_ERROR_CONTEXTS
  #
  # Returns an array of Bio::Gag objects
  #
  # When a block is given, each gag is yielded
  def gags(options={})
    min_disagreeing_proportion = options[:min_disagreeing_proportion]
    min_disagreeing_proportion ||= 0.1
    min_disagreeing_absolute = options[:min_disagreeing_absolute]
    min_disagreeing_absolute ||= 3
    
    options[:acceptable_gag_errors] ||= DEFAULT_GAG_ERROR_CONTEXTS
    
    log = Bio::Log::LoggerPlus['bio-gag']
    
    piles = []
    gags = []
    
    each do |pile|
      options[:progressbar].inc unless options[:progressbar].nil?
      
      if piles.length < 2
        #log.debug "Piles cache for this reference sequence less than length 2"
        piles = [piles, pile].flatten
        next
      elsif piles.length < 3
        #log.debug "Piles cache for this reference sequence becoming full"
        piles = [piles, pile].flatten
      elsif piles[1].ref_name != pile.ref_name
        #log.debug "Piles cache removed - moving to new contig"
        piles = [pile]
        next
      else
        #log.debug "Piles cache regular push through"
        piles = [piles[1], piles[2], pile].flatten
      end
      log.debug "Current piles now at #{piles[0].ref_name}, #{piles.collect{|pile| "#{pile.pos}/#{pile.ref_base}"}.join(', ')}" if log.debug?
      
      # if not at the start/end of the contig
      first = piles[0]
      second = piles[1]
      third = piles[2]
      
      # Require particular sequences in the reference sequence
      ref_bases = "#{first.ref_base.upcase}#{second.ref_base.upcase}#{third.ref_base.upcase}"
      index = options[:acceptable_gag_errors].index(ref_bases)
      if index.nil?
        log.debug "Sequence #{ref_bases} does not match whitelist, so not calling a gag" if log.debug?
        next
      end
      gag_sequence = options[:acceptable_gag_errors][index]
      
      # all reads that have a single insertion after the first or second position, but not both  
      inserting_reads = [first.reads, second.reads].flatten.uniq.select do |read|
        !(read.insertions[first.pos] and read.insertions[second.pos]) and
        (read.insertions[first.pos] or read.insertions[second.pos])
      end
      log.debug "Inserting reads after filtering: #{inserting_reads.inspect}" if log.debug?
      
      # ignore regions that aren't ever going to make it past the next filter
      if inserting_reads.length < min_disagreeing_absolute or inserting_reads.length.to_f/first.coverage < min_disagreeing_proportion
        log.debug "Insufficient disagreement at step 1, so not calling a gag" if log.debug?
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
      log.debug "Direction counts of insertions: #{direction_counts.inspect}" if log.debug?
      log.debug "Base counts of insertions: #{base_counts.inspect}" if log.debug?
      max_direction = direction_counts['+']>direction_counts['-'] ? '+' : '-'
      max_base = base_counts.max do |a,b|
        a[1] <=> b[1]
      end[0]
      log.debug "Picking max direction #{max_direction} and max base #{max_base}" if log.debug?
      
      # Only accept positions that are inserting a single base
      if max_base.length > 1
        log.debug "Maximal insertion is too long, so not calling a gag" if log.debug?
        next
      end
      
      counted_inserts = inserting_reads.select do |read|
        insert = read.insertions[first.pos]
        insert ||= read.insertions[second.pos]
        insert.upcase!
        if read.direction == max_direction and insert == max_base
          # Remove reads that don't match the first and third bases like the consensus sequence
          read.sequence[read.sequence.length-1] == third.ref_base and
          read.sequence[read.sequence.length-3] == first.ref_base
        else
          false
        end
      end
      log.debug "Reads counting after final filtering: #{counted_inserts.inspect}" if log.debug?
      
      coverage = (first.coverage+second.coverage+third.coverage).to_f / 3.0
      coverage_percent = counted_inserts.length.to_f / coverage
      log.debug "Final abundance calculations: max base #{max_base} (comparison base #{second.ref_base.upcase}) occurs #{counted_inserts.length} times compared to coverage #{coverage} (#{coverage_percent*10}%)" if log.debug?
      if max_base != second.ref_base.upcase or # first and second bases must be the same 
        counted_inserts.length < min_disagreeing_absolute or # require 3 bases in that maximal direction
        coverage_percent < min_disagreeing_proportion # at least 10% of reads with disagree with the consensus and agree with the gag
        log.debug "Failed final abundance cutoffs, so not calling a gag" if log.debug?
        next
      end
      
      # alright, gamut navigated. We have a match, record it
      gag = Bio::Gag.new(second.pos, piles, first.ref_name)
      gags.push gag
      log.debug "Yielding gag #{gag.inspect}"
      yield gag if block_given?
    end
    
    return gags
  end
  
  # Given a hash containing sequence identifier => sequences, where both key and value are plain old Ruby strings, return the hash with any GAG errors in the sequences fixed.
  # If the sequence_id_to_gags argument is specified, the gags are not searched from the pileups. If specified, it should be a hash of reference sequence IDs to an array of Bio::Gag objects
  def fix_gags(hash_of_sequence_ids_to_sequence_strings, sequence_id_to_gags={})
    log = Bio::Log::LoggerPlus['bio-gag']
    
    # Get the gags
    if sequence_id_to_gags == {}
      log.info "Predicting gags from the pileup"
      gags do |gag|
        sequence_id_to_gags[gag.ref_name] ||= []
        sequence_id_to_gags[gag.ref_name].push gag
      end
    else
      log.info "Using pre-specified GAG errors"
    end
    log.info "Found #{sequence_id_to_gags.values.flatten.length} gag errors to fix"
    
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
          #log.debug "Attempting to fix gag at position #{gag.position} in sequence #{seq_id}, which is #{seq.length} bases long"
          fixed = fixed+seq[last_gag..(gag.position-1)]
          fixed = fixed+seq[(gag.position-1)..(gag.position-1)]
          last_gag = gag.position
          #log.debug "After fixing gag at position #{gag.position}, fixed sequence is now #{fixed}"
        end
        fixed = fixed+seq[last_gag..(seq.length-1)]
        fixed_sequences[seq_id] = fixed
      end
    end
    
    unless accounted_for_seq_ids.length == sequence_id_to_gags.length
      log.warn "Unexpectedly found GAG errors in sequences that weren't in the sequence that are to be fixed: Found gags in #{sequence_id_to_gags.length}, but only fixed #{accounted_for_seq_ids.length}"
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
  
  # The base to be inserted. May be derived from @gagging_pileups if they have been specified
  attr_writer :inserted_base
  
  def initialize(position, gagging_pileups, ref_name)
    @position = position
    @gagging_pileups = gagging_pileups
    @ref_name = ref_name
  end
  
  # The base to be inserted. May be manually specified in @inserted_base, otherwise it is the ref_base derived from @gagging_pileups at the inserted position
  def inserted_base
    if @inserted_base.nil?
      @gagging_pileups[1].ref_base
    else
      @inserted_base
    end
  end
end


