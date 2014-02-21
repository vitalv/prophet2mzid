require "spect_id_result"
require 'nokogiri'
require 'set'

# A base class for other file formats. Other formats are meant to inherit from this class, thus Format is basically useless by itself.
# Contains some methods that are applicable to all formats.
# Classes that inherit from Format are used as the means of obtaining information from a file to be used in Search2mzIdentML.
#
# @author Jesse Jashinsky (Aug 2010)
class Format
  # @param [String] file the location of the input file
  # @param [String] database the location of the FASTA database that was used by the search engine

  attr_accessor :threshold, :proteins, :peptides, :databaseName, :results, :numberOfSequences

  def initialize(file, database)
    puts "\nPreparing..."
    
    @file = file
    @database = database
    @databaseName = database
    @threshold = 0
    @proteins = []
    @peptides = []
    @results = []    
    @numberOfSequences = 0
    
    check_input #checks whether file and database exist
    
    @missedMappings = Set.new #Set implements a collection of unordered values with no duplicates. hybdrid of Array's intuitive inter-operation facilities and Hash's fast lookup
    
    @non_obo_mappings = Set.new

    @obo = {}
    #oboe.yaml is a custom made mapping (obo_converter.rb creates obo.yaml. Any additional mapping is added to obo.yaml to create oboe.yaml)
    yml = YAML.load_file File.expand_path(File.dirname(__FILE__) + "/../oboe.yaml")
    yml.each {|x| @obo[x[:pepxml_name]] = [x[:id], x[:mzid_name]]}
    #@obo is a Hash with format:
    #["Discoverer MSF", ["MS:1001564", "Discoverer MSF"]]
    # @obo["Discoverer MSF"] =>  ["MS:1001564", "Discoverer MSF"]
  end
  
  
  # Converts calc_neutral_pep_mass to calculatedMassToCharge
  # @param [Float] mass the mass
  # @param [Float] charge the charge
  # @return [Float] the calculatedMassToCharge
  def calMass(mass, charge)
    (mass + (charge.to_f * 1.00727646677)) / charge
  end
  			
  # Converts calc_neutral_pep_mass to experimentalMassToCharge
  # @param [Float] mass the mass
  # @param [Float] charge the charge
  # @param [Float] diff the diff value
  # @return [Float] the experimentalMassToCharge
  def experiMass(mass, charge, diff)
    ((mass + diff) + (charge.to_f * 1.00727646677)) / charge
  end



  def display_missed_mappings
    if !@missedMappings.empty?
      @missedMappings.each do |term|
        puts "WARNING: \"#{term}\" doesn't map to anything in oboe.yaml, and thus won't be displayed in the mzIdentML file."
      end
    end
  end
  


  # Determines the accession number for the given name.
  #
  # @param [String] name the original name of the paramater
  # @return [Aray(String, String)] the number and name
  def findAccession(name)
    if arr = @obo[name]
      arr
    else
#      @missedMappings << name
#      ["", ""]
      @non_obo_mappings << name
      ["", name]
    end
  end
	
  # Conforms score name to mzIdentML format. Will most likely need to be extended.
  #
  # @param [String] name the name of the score
  # @param [String] engine the name of the search engine
  def conformScoreName(name, engine)
    base = 
      case engine
        when "X! Tandem (k-score)"
          "xtandem"
        when "X! Tandem"
          "xtandem"
        when "MASCOT"
          "mascot"
        when "OMSSA"
          "OMSSA"
        when "Tide"
           "sequest"
        when "Phenyx"
          "Phenyx"
        when "SpectraST"
          "SpectraST"
      end
      
    [base, name].join(':')
  end
  
  
  private
  
  def check_input
    raise(ArgumentError, "Invalid input file") if !File.exist?(@file)
    raise(ArgumentError, "Invalid database input") if !File.exist?(@database)
  end
end
