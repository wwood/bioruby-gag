require 'helper'
require 'tempfile'
require 'open3'

class TestBioGag < Test::Unit::TestCase
  should "find_gag" do
    test = "contig00091 4 C 32  ,,..,,......,,,.....,,.,,,,,,,.,  ~~I~~~u~u~t~~~~~~~~~~~~~~~~~~~~~
contig00091 5 G 32  ,,..,,......,,,.....,,.,,,,,,,.,  {{Ii{{iiii@i{{{iiiii{{i{{{{{{{i{
contig00091 6 A 33  ,,.$.+1A,,.+1A.+1A.+1A.+1A.+1A.+1A,,,.+1A.+1A.+1A.+1A.+1A,,.+1A,,,,,,,.+1A,^].  z{D${{$$$$!${{{$$$$${{${{{{{{{${E
contig00091 7 G 32  ,,.,,.....-1G.,,,.....,,.,,,,,,,.,. aaRaaRRRR&RaaaRRRRRaaRaaaaaaaRaU
contig00091 8 G 32  ,,.,,....*.,,,.....,,.,,,,,,,.,.  aaRaaRRRRZRaaaRRRRRaaRaaaaaaaRaa".gsub(/ +/,"\t")
    gags = Bio::DB::PileupIterator.new(test).gags
    assert_equal [6], gags.collect{|g| g.position}
  end

  should "find_gag with first and third bases different, but whitelisted" do
    test = "contig00091 4 C 32  ,,..,,......,,,.....,,.,,,,,,,.,  ~~I~~~u~u~t~~~~~~~~~~~~~~~~~~~~~
contig00091 5 G 32  ,,..,,......,,,.....,,.,,,,,,,.,  {{Ii{{iiii@i{{{iiiii{{i{{{{{{{i{
contig00091 6 C 33  ,,.$.+1C,,.+1C.+1C.+1C.+1C.+1C.+1C,,,.+1C.+1C.+1C.+1C.+1C,,.+1C,,,,,,,.+1C,^].  z{D${{$$$$!${{{$$$$${{${{{{{{{${E
contig00091 7 A 32  ,,.,,.....-1G.,,,.....,,.,,,,,,,.,. aaRaaRRRR&RaaaRRRRRaaRaaaaaaaRaU
contig00091 8 G 32  ,,.,,....*.,,,.....,,.,,,,,,,.,.  aaRaaRRRRZRaaaRRRRRaaRaaaaaaaRaa".gsub(/ +/,"\t")
    gags = Bio::DB::PileupIterator.new(test).gags
    assert_equal [6], gags.collect{|g| g.position}
  end
  
  should "find no gag when XXX" do
    test = "contig00091 4 C 32  ,,..,,......,,,.....,,.,,,,,,,.,  ~~I~~~u~u~t~~~~~~~~~~~~~~~~~~~~~
contig00091 5 G 32  ,,..,,......,,,.....,,.,,,,,,,.,  {{Ii{{iiii@i{{{iiiii{{i{{{{{{{i{
contig00091 6 G 33  ,,.$.+1A,,.+1A.+1A.+1A.+1A.+1A.+1A,,,.+1A.+1A.+1A.+1A.+1A,,.+1A,,,,,,,.+1A,^].  z{D${{$$$$!${{{$$$$${{${{{{{{{${E
contig00091 7 G 32  ,,.,,.....-1G.,,,.....,,.,,,,,,,.,. aaRaaRRRR&RaaaRRRRRaaRaaaaaaaRaU
contig00091 8 G 32  ,,.,,....*.,,,.....,,.,,,,,,,.,.  aaRaaRRRRZRaaaRRRRRaaRaaaaaaaRaa".gsub(/ +/,"\t")
    gags = Bio::DB::PileupIterator.new(test).gags
    assert_equal [], gags.collect{|g| g.position}
  end
  
  should "find no gag with first and third bases are the same but aren't in the whitelist" do
    test = "contig00091 4 C 32  ,,..,,......,,,.....,,.,,,,,,,.,  ~~I~~~u~u~t~~~~~~~~~~~~~~~~~~~~~
contig00091 5 C 32  ,,..,,......,,,.....,,.,,,,,,,.,  {{Ii{{iiii@i{{{iiiii{{i{{{{{{{i{
contig00091 6 A 33  ,,.$.+1A,,.+1A.+1A.+1A.+1A.+1A.+1A,,,.+1A.+1A.+1A.+1A.+1A,,.+1A,,,,,,,.+1A,^].  z{D${{$$$$!${{{$$$$${{${{{{{{{${E
contig00091 7 C 32  ,,.,,.....-1G.,,,.....,,.,,,,,,,.,. aaRaaRRRR&RaaaRRRRRaaRaaaaaaaRaU
contig00091 8 G 32  ,,.,,....*.,,,.....,,.,,,,,,,.,.  aaRaaRRRRZRaaaRRRRRaaRaaaaaaaRaa".gsub(/ +/,"\t")
    gags = Bio::DB::PileupIterator.new(test).gags
    assert_equal [], gags.collect{|g| g.position}
  end
  
  should 'find gags when ref_base is lower-case' do
    test = "contig00036 56618 G 19  ......-1C,...,.....,,.  ~~~~~il~~~l~~~~~ll~
contig00036 56619 c 19  .....*,+1c...,+1c.....,+1c,+1c. tulutI$uuu$uuuuu$$u
contig00036 56620 T 21  ......,...,.....,,.^].^], qrirqILrrrLrrrrrLLrNE".gsub(/ +/,"\t")
    
    gags = Bio::DB::PileupIterator.new(test).gags
    assert_equal [56619], gags.collect{|g| g.position}
  end

  should "fix gag" do
    test = "contig00091 1 G 32  ,,..,,......,,,.....,,.,,,,,,,.,  {;c{{{l{l{l{{{{{{{{{{{{{{{{{{{{U
contig00091 2 T 32  ,,.-1T.,,.-1T..-1T..-1T.,,,.....,,.,,,,,,,.,  a`$aaa!a!a!aaaaaaaaaaaaaaaaaaaaa
contig00091 3 T 32  ,,*.,,*.*.*.,,,.....,,.,,,,,,,.,  a`Iaaauauataaaaaaaaaaaaaaaaaaaaa
contig00091 4 C 32  ,,..,,......,,,.....,,.,,,,,,,.,  ~~I~~~u~u~t~~~~~~~~~~~~~~~~~~~~~
contig00091 5 G 32  ,,..,,......,,,.....,,.,,,,,,,.,  {{Ii{{iiii@i{{{iiiii{{i{{{{{{{i{
contig00091 6 A 33  ,,.$.+1A,,.+1A.+1A.+1A.+1A.+1A.+1A,,,.+1A.+1A.+1A.+1A.+1A,,.+1A,,,,,,,.+1A,^].  z{D${{$$$$!${{{$$$$${{${{{{{{{${E
contig00091 7 G 32  ,,.,,.....-1G.,,,.....,,.,,,,,,,.,. aaRaaRRRR&RaaaRRRRRaaRaaaaaaaRaU
contig00091 8 G 32  ,,.,,....*.,,,.....,,.,,,,,,,.,.  aaRaaRRRRZRaaaRRRRRaaRaaaaaaaRaa
contig00091 9 C 32  ,,.,,......,,,.....,,.,,,,,,,.,.  ~~i~~~~~~Z~~~~~~~~~~~~~~~~~~~~~r
contig00091 10  A 33  ,,.,,......,,,.....,,.,,,,,,,.,.^]. aaPaa^aaaYaaaaaaaaaaaaaaaaaaaaaaB".gsub(/ +/,"\t")
    hash = {'contig00091' => 'GTTCGAGGC'}
    expe = {'contig00091' => 'GTTCGAAGGC'}
    assert_equal expe, Bio::DB::PileupIterator.new(test).fix_gags(hash)
  end
  
  should 'fix a gag when ref_base is lower-case' do
    test = "contig00036 1 G 19  ......-1C,...,.....,,.  ~~~~~il~~~l~~~~~ll~
contig00036 2 c 19  .....*,+1c...,+1c.....,+1c,+1c. tulutI$uuu$uuuuu$$u
contig00036 3 T 21  ......,...,.....,,.^].^], qrirqILrrrLrrrrrLLrNE".gsub(/ +/,"\t")
    
    hash = {'contig00036' => 'GCT'}
    expe = {'contig00036' => 'GCCT'}
    assert_equal expe, Bio::DB::PileupIterator.new(test).fix_gags(hash)
  end
  
  should "fix gag prespecified" do
    test = "contig00091 1 G 32  ,,..,,......,,,.....,,.,,,,,,,.,  {;c{{{l{l{l{{{{{{{{{{{{{{{{{{{{U
contig00091 2 T 32  ,,.-1T.,,.-1T..-1T..-1T.,,,.....,,.,,,,,,,.,  a`$aaa!a!a!aaaaaaaaaaaaaaaaaaaaa
contig00091 3 T 32  ,,*.,,*.*.*.,,,.....,,.,,,,,,,.,  a`Iaaauauataaaaaaaaaaaaaaaaaaaaa
contig00091 4 C 32  ,,..,,......,,,.....,,.,,,,,,,.,  ~~I~~~u~u~t~~~~~~~~~~~~~~~~~~~~~
contig00091 5 G 32  ,,..,,......,,,.....,,.,,,,,,,.,  {{Ii{{iiii@i{{{iiiii{{i{{{{{{{i{
contig00091 6 A 33  ,,.$.+1A,,.+1A.+1A.+1A.+1A.+1A.+1A,,,.+1A.+1A.+1A.+1A.+1A,,.+1A,,,,,,,.+1A,^].  z{D${{$$$$!${{{$$$$${{${{{{{{{${E
contig00091 7 G 32  ,,.,,.....-1G.,,,.....,,.,,,,,,,.,. aaRaaRRRR&RaaaRRRRRaaRaaaaaaaRaU
contig00091 8 G 32  ,,.,,....*.,,,.....,,.,,,,,,,.,.  aaRaaRRRRZRaaaRRRRRaaRaaaaaaaRaa
contig00091 9 C 32  ,,.,,......,,,.....,,.,,,,,,,.,.  ~~i~~~~~~Z~~~~~~~~~~~~~~~~~~~~~r
contig00091 10  A 33  ,,.,,......,,,.....,,.,,,,,,,.,.^]. aaPaa^aaaYaaaaaaaaaaaaaaaaaaaaaaB".gsub(/ +/,"\t")
    hash = {'contig00091' => 'GTTCGAGGC'}
    expe = {'contig00091' => 'GTTTCGAGGC'}
    gag1 = Bio::Gag.new(2,nil,'contig00091')
    gags = {'contig00091' => [gag1]}
    assert_equal expe, gags = Bio::DB::PileupIterator.new(test).fix_gags(hash, gags)
  end
  
  should "fix gag prespecified in 2 seqs" do
    hash = {'contig00091' => 'GTTCGAGGC',
      'contig00092' => 'GAGTTCGAGGC'}
    expe = {'contig00091' => 'GTTTCGAGGC',
      'contig00092' => 'GAGTTCGAGGC'}
      
    gag1 = Bio::Gag.new(2,nil,'contig00091')
    gags = {'contig00091' => [gag1]}
    assert_equal expe, gags = Bio::DB::PileupIterator.new('').fix_gags(hash, gags)
    
    gag2 = Bio::Gag.new(8,nil,'contig00092')
    gags = {'contig00091' => [gag1], 'contig00092' => [gag2]}
    expe = {'contig00091' => 'GTTTCGAGGC',
      'contig00092' => 'GAGTTCGAAGGC'}
    assert_equal expe, gags = Bio::DB::PileupIterator.new('').fix_gags(hash, gags)
  end

  should "fix 2 gags" do
    test = "contig00091 1 G 32  ,,..,,......,,,.....,,.,,,,,,,.,  {;c{{{l{l{l{{{{{{{{{{{{{{{{{{{{U
contig00091 2 T 32  ,,.-1T.,,.-1T..-1T..-1T.,,,.....,,.,,,,,,,.,  a`$aaa!a!a!aaaaaaaaaaaaaaaaaaaaa
contig00091 3 T 32  ,,*.,,*.*.*.,,,.....,,.,,,,,,,.,  a`Iaaauauataaaaaaaaaaaaaaaaaaaaa
contig00091 4 C 32  ,,..,,......,,,.....,,.,,,,,,,.,  ~~I~~~u~u~t~~~~~~~~~~~~~~~~~~~~~
contig00091 5 G 32  ,,..,,......,,,.....,,.,,,,,,,.,  {{Ii{{iiii@i{{{iiiii{{i{{{{{{{i{
contig00091 6 A 33  ,,..+1A,,.+1A.+1A.+1A.+1A.+1A.+1A,,,.+1A.+1A.+1A.+1A.+1A,,.+1A,,,,,,,.+1A,^].  z{D${{$$$$!${{{$$$$${{${{{{{{{${E
contig00091 7 G 32  ,,..,,.....-1G.,,,.....,,.,,,,,,,.,. aaRaaRRRR&RaaaRRRRRaaRaaaaaaaRaU
contig00091 8 G 32  ,,..,,......,,,.....,,.,,,,,,,.,.  {{Ii{{iiii@i{{{iiiii{{i{{{{{{{i{
contig00091 9 A 33  ,,.$.+1A,,.+1A.+1A.+1A.+1A.+1A.+1A,,,.+1A.+1A.+1A.+1A.+1A,,.+1A,,,,,,,.+1A,.  z{D${{$$$$!${{{$$$$${{${{{{{{{${E
contig00091 10 G 32  ,,.,,.....-1G.,,,.....,,.,,,,,,,.,. aaRaaRRRR&RaaaRRRRRaaRaaaaaaaRaU
contig00091 11 G 32  ,,.,,....*.,,,.....,,.,,,,,,,.,.  aaRaaRRRRZRaaaRRRRRaaRaaaaaaaRaa
contig00091 12 C 32  ,,.,,......,,,.....,,.,,,,,,,.,.  ~~i~~~~~~Z~~~~~~~~~~~~~~~~~~~~~r
contig00091 13  A 33  ,,.,,......,,,.....,,.,,,,,,,.,.^]. aaPaa^aaaYaaaaaaaaaaaaaaaaaaaaaaB".gsub(/ +/,"\t")

    hash = {'contig00091' => 'GTTCGAGGAGGCA'}
    expe = {'contig00091' => 'GTTCGAAGGAAGGCA'}
    assert_equal expe, gags = Bio::DB::PileupIterator.new(test).fix_gags(hash)
  end

  should "run gagger predict ok" do
    test = "contig00091 4 C 32  ,,..,,......,,,.....,,.,,,,,,,.,  ~~I~~~u~u~t~~~~~~~~~~~~~~~~~~~~~
contig00091 5 G 32  ,,..,,......,,,.....,,.,,,,,,,.,  {{Ii{{iiii@i{{{iiiii{{i{{{{{{{i{
contig00091 6 A 33  ,,.$.+1A,,.+1A.+1A.+1A.+1A.+1A.+1A,,,.+1A.+1A.+1A.+1A.+1A,,.+1A,,,,,,,.+1A,^].  z{D${{$$$$!${{{$$$$${{${{{{{{{${E
contig00091 7 G 32  ,,.,,.....-1G.,,,.....,,.,,,,,,,.,. aaRaaRRRR&RaaaRRRRRaaRaaaaaaaRaU
contig00091 8 G 32  ,,.,,....*.,,,.....,,.,,,,,,,.,.  aaRaaRRRRZRaaaRRRRRaaRaaaaaaaRaa".gsub(/ +/,"\t")

    Tempfile.open('stdin') do |input|
      input.puts test
      input.close
      
      command = File.join([File.dirname(__FILE__),%w(.. bin gag)].flatten)+ ' --trace warn --no-progress'+' -p '+input.path
      out = nil
      err = nil
      Open3.popen3(command) do |stdin, stdout, stderr|
        out = stdout.readlines
        err = stderr.readlines
      end
      assert_equal [], err
      assert_equal [
        "ref_name\tposition\tinserted_base\tcontext\n",
        "contig00091\t6\tA\tGAG\n"
        ], out
    end
  end
  
  should 'run gagger predict with specified contexts correctly' do
    test = "contig00091 4 C 32  ,,..,,......,,,.....,,.,,,,,,,.,  ~~I~~~u~u~t~~~~~~~~~~~~~~~~~~~~~
contig00091 5 G 32  ,,..,,......,,,.....,,.,,,,,,,.,  {{Ii{{iiii@i{{{iiiii{{i{{{{{{{i{
contig00091 6 A 33  ,,.$.+1A,,.+1A.+1A.+1A.+1A.+1A.+1A,,,.+1A.+1A.+1A.+1A.+1A,,.+1A,,,,,,,.+1A,^].  z{D${{$$$$!${{{$$$$${{${{{{{{{${E
contig00091 7 G 32  ,,.,,.....-1G.,,,.....,,.,,,,,,,.,. aaRaaRRRR&RaaaRRRRRaaRaaaaaaaRaU
contig00091 8 G 32  ,,.,,....*.,,,.....,,.,,,,,,,.,.  aaRaaRRRRZRaaaRRRRRaaRaaaaaaaRaa".gsub(/ +/,"\t")

    Tempfile.open('stdin') do |input|
      input.puts test
      input.close
      
      command = File.join([File.dirname(__FILE__),%w(.. bin gag)].flatten)+ ' --trace warn --contexts GAG,CGC --no-progress'+' -p '+input.path
      out = nil
      err = nil
      Open3.popen3(command) do |stdin, stdout, stderr|
        out = stdout.readlines
        err = stderr.readlines
      end
      assert_equal [], err
      assert_equal [
        "ref_name\tposition\tinserted_base\tcontext\n",
        "contig00091\t6\tA\tGAG\n"
        ], out
        
        
      command = File.join([File.dirname(__FILE__),%w(.. bin gag)].flatten)+ ' --trace warn --contexts TGC,GCA --no-progress'+' -p '+input.path
      out = nil
      err = nil
      Open3.popen3(command) do |stdin, stdout, stderr|
        out = stdout.readlines
        err = stderr.readlines
      end
      assert_equal [], err
      assert_equal [
        "ref_name\tposition\tinserted_base\tcontext\n",
        # "contig00091\t6\tA\tGAG\n"
        ], out
    end
  end

  should "run gagger fix ok without gags pre-specified" do
    test = "contig00091 1 G 32  ,,..,,......,,,.....,,.,,,,,,,.,  {;c{{{l{l{l{{{{{{{{{{{{{{{{{{{{U
contig00091 2 T 32  ,,.-1T.,,.-1T..-1T..-1T.,,,.....,,.,,,,,,,.,  a`$aaa!a!a!aaaaaaaaaaaaaaaaaaaaa
contig00091 3 T 32  ,,*.,,*.*.*.,,,.....,,.,,,,,,,.,  a`Iaaauauataaaaaaaaaaaaaaaaaaaaa
contig00091 4 C 32  ,,..,,......,,,.....,,.,,,,,,,.,  ~~I~~~u~u~t~~~~~~~~~~~~~~~~~~~~~
contig00091 5 G 32  ,,..,,......,,,.....,,.,,,,,,,.,  {{Ii{{iiii@i{{{iiiii{{i{{{{{{{i{
contig00091 6 A 33  ,,..+1A,,.+1A.+1A.+1A.+1A.+1A.+1A,,,.+1A.+1A.+1A.+1A.+1A,,.+1A,,,,,,,.+1A,^].  z{D${{$$$$!${{{$$$$${{${{{{{{{${E
contig00091 7 G 32  ,,..,,.....-1G.,,,.....,,.,,,,,,,.,. aaRaaRRRR&RaaaRRRRRaaRaaaaaaaRaU
contig00091 8 G 32  ,,..,,......,,,.....,,.,,,,,,,.,.  {{Ii{{iiii@i{{{iiiii{{i{{{{{{{i{
contig00091 9 A 33  ,,.$.+1A,,.+1A.+1A.+1A.+1A.+1A.+1A,,,.+1A.+1A.+1A.+1A.+1A,,.+1A,,,,,,,.+1A,.  z{D${{$$$$!${{{$$$$${{${{{{{{{${E
contig00091 10 G 32  ,,.,,.....-1G.,,,.....,,.,,,,,,,.,. aaRaaRRRR&RaaaRRRRRaaRaaaaaaaRaU
contig00091 11 G 32  ,,.,,....*.,,,.....,,.,,,,,,,.,.  aaRaaRRRRZRaaaRRRRRaaRaaaaaaaRaa
contig00091 12 C 32  ,,.,,......,,,.....,,.,,,,,,,.,.  ~~i~~~~~~Z~~~~~~~~~~~~~~~~~~~~~r
contig00091 13  A 33  ,,.,,......,,,.....,,.,,,,,,,.,.^]. aaPaa^aaaYaaaaaaaaaaaaaaaaaaaaaaB".gsub(/ +/,"\t")
    Tempfile.open('test_gag_fix') do |tempfile|
      tempfile.puts '>contig00091'
      tempfile.puts 'GTTCGAGGAGGCA'
      tempfile.close
      
      Tempfile.open('stdin') do |input|
        input.puts test
        input.close

        command = File.join([File.dirname(__FILE__),%w(.. bin gag)].flatten)+' --trace error --fix '+tempfile.path+' -p '+input.path
        out = nil
        err = nil
        Open3.popen3(command) do |stdin, stdout, stderr|
          out = stdout.readlines
          err = stderr.readlines
        end
        assert_equal [], err
        assert_equal [
        ">contig00091\n",
        "GTTCGAAGGAAGGCA\n"
        ], out
      end
    end
  end

  should "run gagger fix ok with fasta comments" do
    test = "contig00091 1 G 32  ,,..,,......,,,.....,,.,,,,,,,.,  {;c{{{l{l{l{{{{{{{{{{{{{{{{{{{{U
contig00091 2 T 32  ,,.-1T.,,.-1T..-1T..-1T.,,,.....,,.,,,,,,,.,  a`$aaa!a!a!aaaaaaaaaaaaaaaaaaaaa
contig00091 3 T 32  ,,*.,,*.*.*.,,,.....,,.,,,,,,,.,  a`Iaaauauataaaaaaaaaaaaaaaaaaaaa
contig00091 4 C 32  ,,..,,......,,,.....,,.,,,,,,,.,  ~~I~~~u~u~t~~~~~~~~~~~~~~~~~~~~~
contig00091 5 G 32  ,,..,,......,,,.....,,.,,,,,,,.,  {{Ii{{iiii@i{{{iiiii{{i{{{{{{{i{
contig00091 6 A 33  ,,..+1A,,.+1A.+1A.+1A.+1A.+1A.+1A,,,.+1A.+1A.+1A.+1A.+1A,,.+1A,,,,,,,.+1A,^].  z{D${{$$$$!${{{$$$$${{${{{{{{{${E
contig00091 7 G 32  ,,..,,.....-1G.,,,.....,,.,,,,,,,.,. aaRaaRRRR&RaaaRRRRRaaRaaaaaaaRaU
contig00091 8 G 32  ,,..,,......,,,.....,,.,,,,,,,.,.  {{Ii{{iiii@i{{{iiiii{{i{{{{{{{i{
contig00091 9 A 33  ,,.$.+1A,,.+1A.+1A.+1A.+1A.+1A.+1A,,,.+1A.+1A.+1A.+1A.+1A,,.+1A,,,,,,,.+1A,.  z{D${{$$$$!${{{$$$$${{${{{{{{{${E
contig00091 10 G 32  ,,.,,.....-1G.,,,.....,,.,,,,,,,.,. aaRaaRRRR&RaaaRRRRRaaRaaaaaaaRaU
contig00091 11 G 32  ,,.,,....*.,,,.....,,.,,,,,,,.,.  aaRaaRRRRZRaaaRRRRRaaRaaaaaaaRaa
contig00091 12 C 32  ,,.,,......,,,.....,,.,,,,,,,.,.  ~~i~~~~~~Z~~~~~~~~~~~~~~~~~~~~~r
contig00091 13  A 33  ,,.,,......,,,.....,,.,,,,,,,.,.^]. aaPaa^aaaYaaaaaaaaaaaaaaaaaaaaaaB".gsub(/ +/,"\t")
    Tempfile.open('test_gag_fix') do |tempfile|
      tempfile.puts '>contig00091 with comment'
      tempfile.puts 'GTTCGAGGAGGCA'
      tempfile.close

      Tempfile.open('stdin') do |input|
        input.puts test
        input.close
        command = File.join([File.dirname(__FILE__),%w(.. bin gag)].flatten)+' --trace error --fix '+tempfile.path+' -p '+input.path
      
        out = nil
        err = nil
        Open3.popen3(command) do |stdin, stdout, stderr|
          out = stdout.readlines
          err = stderr.readlines
        end
        assert_equal [], err
        assert_equal [
        ">contig00091\n",
        "GTTCGAAGGAAGGCA\n"
        ], out
      end
    end
  end

  should "run gagger fix when some sequences don't have gag errors" do
    test = "contig00091 1 C 32  ,,..,,......,,,.....,,.,,,,,,,.,  ~~I~~~u~u~t~~~~~~~~~~~~~~~~~~~~~
contig00091 2 G 32  ,,..,,......,,,.....,,.,,,,,,,.,  {{Ii{{iiii@i{{{iiiii{{i{{{{{{{i{
contig00091 3 A 33  ,,.$.+1A,,.+1A.+1A.+1A.+1A.+1A.+1A,,,.+1A.+1A.+1A.+1A.+1A,,.+1A,,,,,,,.+1A,^].  z{D${{$$$$!${{{$$$$${{${{{{{{{${E
contig00091 4 G 32  ,,.,,.....-1G.,,,.....,,.,,,,,,,.,. aaRaaRRRR&RaaaRRRRRaaRaaaaaaaRaU
contig00091 5 G 32  ,,.,,....*.,,,.....,,.,,,,,,,.,.  aaRaaRRRRZRaaaRRRRRaaRaaaaaaaRaa".gsub(/ +/,"\t")

    Tempfile.open('test_gag_fix') do |tempfile|
      tempfile.puts '>contig00091 with comment'
      tempfile.puts 'CGAGG'
      tempfile.puts '>contig00092'
      tempfile.puts 'ATGC'
      tempfile.close

      Tempfile.open('stdin') do |input|
        input.puts test
        input.close
        command = File.join([File.dirname(__FILE__),%w(.. bin gag)].flatten)+' --trace error --fix '+tempfile.path+' -p '+input.path
        out = nil
        err = nil
        Open3.popen3(command) do |stdin, stdout, stderr|
          out = stdout.readlines
          err = stderr.readlines
        end
        assert_equal [], err
        assert_equal [
        ">contig00091\n",
        "CGAAGG\n",
        ">contig00092\n",
        "ATGC\n"
        ], out
      end
    end
  end


  should "run gagger fix ok, but warn, when there's less sequences than gag errors" do
    test = "contig00091 1 C 32  ,,..,,......,,,.....,,.,,,,,,,.,  ~~I~~~u~u~t~~~~~~~~~~~~~~~~~~~~~
contig00091 2 G 32  ,,..,,......,,,.....,,.,,,,,,,.,  {{Ii{{iiii@i{{{iiiii{{i{{{{{{{i{
contig00091 3 A 33  ,,.$.+1A,,.+1A.+1A.+1A.+1A.+1A.+1A,,,.+1A.+1A.+1A.+1A.+1A,,.+1A,,,,,,,.+1A,^].  z{D${{$$$$!${{{$$$$${{${{{{{{{${E
contig00091 4 G 32  ,,.,,.....-1G.,,,.....,,.,,,,,,,.,. aaRaaRRRR&RaaaRRRRRaaRaaaaaaaRaU
contig00091 5 G 32  ,$,$.$,$,$.$.$.$.$*$.$,$,$,$.$.$.$.$.$,$,$.$,$,$,$,$,$,$,$.$,$.$  aaRaaRRRRZRaaaRRRRRaaRaaaaaaaRaa
contig00090 1 C 32  ,,..,,......,,,.....,,.,,,,,,,.,  ~~I~~~u~u~t~~~~~~~~~~~~~~~~~~~~~
contig00090 2 G 32  ,,..,,......,,,.....,,.,,,,,,,.,  {{Ii{{iiii@i{{{iiiii{{i{{{{{{{i{
contig00090 3 A 33  ,,.$.+1A,,.+1A.+1A.+1A.+1A.+1A.+1A,,,.+1A.+1A.+1A.+1A.+1A,,.+1A,,,,,,,.+1A,^].  z{D${{$$$$!${{{$$$$${{${{{{{{{${E
contig00090 4 G 32  ,,.,,.....-1G.,,,.....,,.,,,,,,,.,. aaRaaRRRR&RaaaRRRRRaaRaaaaaaaRaU
contig00090 5 G 32  ,,.,,....*.,,,.....,,.,,,,,,,.,.  aaRaaRRRRZRaaaRRRRRaaRaaaaaaaRaa".gsub(/ +/,"\t")

    Tempfile.open('test_gag_fix') do |tempfile|
      tempfile.puts '>contig00091 with comment'
      tempfile.puts 'CGAGG'
      tempfile.close

      Tempfile.open('stdin') do |input|
        input.puts test
        input.close
        
        command = File.join([File.dirname(__FILE__),%w(.. bin gag)].flatten)+ ' --trace warn --fix '+tempfile.path+' -p '+input.path
        out = nil
        err = nil
        Open3.popen3(command) do |stdin, stdout, stderr|
          stdin.puts test
          stdin.close
          out = stdout.readlines
          err = stderr.readlines
        end
        assert_equal [" WARN bio-gag: Unexpectedly found GAG errors in sequences that weren't in the sequence that are to be fixed: Found gags in 2, but only fixed 1\n"], err
        assert_equal [
        ">contig00091\n",
        "CGAAGG\n",
        ], out
      end
    end
  end
  
  should "run gagger with --debug without any big problems" do
    test = "contig00091 1 C 32  ,,..,,......,,,.....,,.,,,,,,,.,  ~~I~~~u~u~t~~~~~~~~~~~~~~~~~~~~~
contig00091 2 G 32  ,,..,,......,,,.....,,.,,,,,,,.,  {{Ii{{iiii@i{{{iiiii{{i{{{{{{{i{
contig00091 3 A 33  ,,.$.+1A,,.+1A.+1A.+1A.+1A.+1A.+1A,,,.+1A.+1A.+1A.+1A.+1A,,.+1A,,,,,,,.+1A,^].  z{D${{$$$$!${{{$$$$${{${{{{{{{${E
contig00091 4 G 32  ,,.,,.....-1G.,,,.....,,.,,,,,,,.,. aaRaaRRRR&RaaaRRRRRaaRaaaaaaaRaU
contig00091 5 G 32  ,$,$.$,$,$.$.$.$.$*$.$,$,$,$.$.$.$.$.$,$,$.$,$,$,$,$,$,$,$.$,$.$  aaRaaRRRRZRaaaRRRRRaaRaaaaaaaRaa
contig00090 1 C 32  ,,..,,......,,,.....,,.,,,,,,,.,  ~~I~~~u~u~t~~~~~~~~~~~~~~~~~~~~~
contig00090 2 G 32  ,,..,,......,,,.....,,.,,,,,,,.,  {{Ii{{iiii@i{{{iiiii{{i{{{{{{{i{
contig00090 3 A 33  ,,.$.+1A,,.+1A.+1A.+1A.+1A.+1A.+1A,,,.+1A.+1A.+1A.+1A.+1A,,.+1A,,,,,,,.+1A,^].  z{D${{$$$$!${{{$$$$${{${{{{{{{${E
contig00090 4 G 32  ,,.,,.....-1G.,,,.....,,.,,,,,,,.,. aaRaaRRRR&RaaaRRRRRaaRaaaaaaaRaU
contig00090 5 G 32  ,,.,,....*.,,,.....,,.,,,,,,,.,.  aaRaaRRRRZRaaaRRRRRaaRaaaaaaaRaa".gsub(/ +/,"\t")

    Tempfile.open('test_gag_fix') do |tempfile|
      tempfile.puts '>contig00091 with comment'
      tempfile.puts 'CGAGG'
      tempfile.close

      Tempfile.open('stdin') do |input|
        input.puts test
        input.close
        
        command = File.join([File.dirname(__FILE__),%w(.. bin gag)].flatten)+ ' --trace debug --fix '+tempfile.path+' -p '+input.path
        out = nil
        err = nil
        Open3.popen3(command) do |stdin, stdout, stderr|
          out = stdout.readlines
          err = stderr.readlines
        end
        assert err.length > 1, "expected more errors"
        assert_equal [
        ">contig00091\n",
        "CGAAGG\n",
        ], out
      end
    end
  end

 should "run gagger fix ok with prespecified gags" do
    test = ""
    Tempfile.open('test_gag_fix') do |tempfile|
      tempfile.puts '>contig00091'
      tempfile.puts 'GTTCGAGGAGGCA'
      tempfile.close

      Tempfile.open('gags_prespecified') do |gags_file|
        gags_file.puts %w(ref_name  position  inserted_base context).join("\t")
        gags_file.puts %w(contig00091 2 G CTC).join("\t")
        gags_file.puts %w(contig00091 4 G CTC).join("\t")
        gags_file.close

        Tempfile.open('stdin') do |input|
          input.puts test
          input.close
          
          command = File.join([File.dirname(__FILE__),%w(.. bin gag)].flatten)+" --trace error --fix #{tempfile.path} --gags #{gags_file.path} -p "+input.path
          out = nil
          err = nil
          Open3.popen3(command) do |stdin, stdout, stderr|
            out = stdout.readlines
            err = stderr.readlines
          end
          assert_equal [], err
          assert_equal [
            ">contig00091\n",
            "GTTTCCGAGGAGGCA\n"
            ], out
        end
      end
    end
  end

end
