#!/usr/bin/ruby

require "rubygems"
require "bundler/setup"
# require 'optparse'

require "#{File.expand_path(File.dirname(__FILE__))}/../lib/search2mzidentml.rb" 
#__FILE__ is a reference to current file (here it's search2mzidentml_cl.rb


if ARGV.size != 3 #2 
  puts "\nusage: #{File.basename(__FILE__)} pepxmlFile protxmlFile database"
  puts "pepxmlFile: The location of the pep.xml file "
  puts "protxmlFile: The location of the prot.xml file "
  puts "database: The location of the FASTA database\n\n"
  exit
end

#Vital:
#In order to combine data from  pepxml and protxml into mzIdentML, 3 ARGV are needed now
# if ARGV.size != 3 #check order of ARGV too!
#   puts "\nusage: #{File.basename(__FILE__)} inputFile.pep.xml inputFILE.prot.xml database"
#   exit
# end


begin

  pepxml = PepXML.new(ARGV[0], ARGV[2]) 
  #pepxml = PepXML.new(ARGV[0], ARGV[1])
  protxml = ProtXML.new(ARGV[1], ARGV[2])

  Search2mzIdentML.new(pepxml, protxml).convert
  #Search2mzIdentML.new(pepxml).convert

rescue Exception => msg
#  $stderr.print "\n\tError: #{$!}\n"
	#$! last error message
  puts msg
  puts msg.backtrace
  puts msg.inspect
end
