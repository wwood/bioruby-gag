
require 'bio-logger'
Bio::Log::LoggerPlus.new('bio-gag')

$:.unshift File.join(File.dirname(__FILE__),'../../bioruby-pileup_iterator_no_debug/lib/')
require 'bio-pileup_iterator'
require 'bio/db/gag'

