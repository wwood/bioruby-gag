require 'helper'
require 'tempfile'
require 'open3'

class TestBioGag < Test::Unit::TestCase
  should "find_gag" do
    test = "contig00091 4 C 32  ,,..,,......,,,.....,,.,,,,,,,.,  ~~I~~~u~u~t~~~~~~~~~~~~~~~~~~~~~
contig00091 5 G 32  ,,..,,......,,,.....,,.,,,,,,,.,  {{Ii{{iiii@i{{{iiiii{{i{{{{{{{i{
contig00091 6 A 33  ,,.$.+1A,,.+1A.+1A.+1A.+1A.+1A.+1A,,,.+1A.+1A.+1A.+1A.+1A,,.+1A,,,,,,,.+1A,^].  z{D${{$$$$!${{{$$$$${{${{{{{{{${E
contig00091 7 G 32  ,,.,,.....-1G.,,,.....,,.,,,,,,,.,. aaRaaRRRR&RaaaRRRRRaaRaaaaaaaRaU
contig00091 8 G 32  ,,.,,....*.,,,.....,,.,,,,,,,.,.  aaRaaRRRRZRaaaRRRRRaaRaaaaaaaRaa"
    gags = Bio::DB::PileupIterator.new(test).gags
    assert_equal [6], gags.collect{|g| g.position}
  end

  should "find no gag when XXX" do
    test = "contig00091 4 C 32  ,,..,,......,,,.....,,.,,,,,,,.,  ~~I~~~u~u~t~~~~~~~~~~~~~~~~~~~~~
contig00091 5 G 32  ,,..,,......,,,.....,,.,,,,,,,.,  {{Ii{{iiii@i{{{iiiii{{i{{{{{{{i{
contig00091 6 G 33  ,,.$.+1A,,.+1A.+1A.+1A.+1A.+1A.+1A,,,.+1A.+1A.+1A.+1A.+1A,,.+1A,,,,,,,.+1A,^].  z{D${{$$$$!${{{$$$$${{${{{{{{{${E
contig00091 7 G 32  ,,.,,.....-1G.,,,.....,,.,,,,,,,.,. aaRaaRRRR&RaaaRRRRRaaRaaaaaaaRaU
contig00091 8 G 32  ,,.,,....*.,,,.....,,.,,,,,,,.,.  aaRaaRRRRZRaaaRRRRRaaRaaaaaaaRaa"
    gags = Bio::DB::PileupIterator.new(test).gags
    assert_equal [], gags.collect{|g| g.position}
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
contig00091 10  A 33  ,,.,,......,,,.....,,.,,,,,,,.,.^]. aaPaa^aaaYaaaaaaaaaaaaaaaaaaaaaaB"
    hash = {'contig00091' => 'GTTCGAGGC'}
    expe = {'contig00091' => 'GTTCGAAGGC'}
    assert_equal expe, gags = Bio::DB::PileupIterator.new(test).fix_gags(hash)
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
contig00091 13  A 33  ,,.,,......,,,.....,,.,,,,,,,.,.^]. aaPaa^aaaYaaaaaaaaaaaaaaaaaaaaaaB"

    hash = {'contig00091' => 'GTTCGAGGAGGCA'}
    expe = {'contig00091' => 'GTTCGAAGGAAGGCA'}
    assert_equal expe, gags = Bio::DB::PileupIterator.new(test).fix_gags(hash)
  end

  should "run gagger predict ok" do
    test = "contig00091 4 C 32  ,,..,,......,,,.....,,.,,,,,,,.,  ~~I~~~u~u~t~~~~~~~~~~~~~~~~~~~~~
contig00091 5 G 32  ,,..,,......,,,.....,,.,,,,,,,.,  {{Ii{{iiii@i{{{iiiii{{i{{{{{{{i{
contig00091 6 A 33  ,,.$.+1A,,.+1A.+1A.+1A.+1A.+1A.+1A,,,.+1A.+1A.+1A.+1A.+1A,,.+1A,,,,,,,.+1A,^].  z{D${{$$$$!${{{$$$$${{${{{{{{{${E
contig00091 7 G 32  ,,.,,.....-1G.,,,.....,,.,,,,,,,.,. aaRaaRRRR&RaaaRRRRRaaRaaaaaaaRaU
contig00091 8 G 32  ,,.,,....*.,,,.....,,.,,,,,,,.,.  aaRaaRRRRZRaaaRRRRRaaRaaaaaaaRaa"
    command = File.join([File.dirname(__FILE__),%w(.. bin gag)].flatten)+ ' --debug'
    out = nil
    err = nil
    Open3.popen3(command) do |stdin, stdout, stderr|
      stdin.puts test
      stdin.close
      out = stdout.readlines
      err = stderr.readlines
    end
    assert_equal [], err
    assert_equal [
      "ref_name\tposition\tinserted_base\tcontext\n",
      "contig00091\t6\tA\tGAG\n"
      ], out
  end

  should "run gagger fix ok" do
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
contig00091 13  A 33  ,,.,,......,,,.....,,.,,,,,,,.,.^]. aaPaa^aaaYaaaaaaaaaaaaaaaaaaaaaaB"
    Tempfile.open('test_gag_fix') do |tempfile|
      tempfile.puts '>contig00091'
      tempfile.puts 'GTTCGAGGAGGCA'
      tempfile.close
      
      command = File.join([File.dirname(__FILE__),%w(.. bin gag)].flatten)+' --debug --fix '+tempfile.path
      out = nil
      err = nil
      Open3.popen3(command) do |stdin, stdout, stderr|
        stdin.puts test
        stdin.close
        out = stdout.readlines
        err = stderr.readlines
      end
      assert_equal [], err
      assert_equal [
      "'>contig00091'\n",
      "GTTCGAAGGAAGGCA\n"
      ], out
    end
  end
end
