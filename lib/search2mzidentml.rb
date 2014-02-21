#$: contains the directories that make up the path Ruby searches when loading an external file
#so here I'm including all files in .
$: << "#{File.expand_path(File.dirname(__FILE__))}/"

require "formats/pepxml"
require "formats/protxml"
require 'nokogiri'

# Creates an mzIdentML file from a file type created by a search engine, using the format classes such as PepXML.
#
# @author Jesse Jashinsky (Aug 2010)
#
# Added functionality: Uses data from pepXML and protXML combining both to create mzIdentML (mzIdentML1_1_0 element <ProteinDetectionList> is now included)
# Author: Vital Vialas (Apr 2012)

@@beginning_time = Time.now
(1..10000).each { |i| i}

#~ puts "Time elapsed: #{end_time - beginning_time} seconds"


class Search2mzIdentML
  # @param [String] format a Format object

  #def initialize(pepxml)
  def initialize(pepxml, protxml)
    @pepxml = pepxml
    @protxml = protxml
  end
  
  # Starts the Nokogiri build process. Other methods build the different parts of the file. Root is depth 0
  #
  # @option opts [] currently does nothing
  def convert(opts={})
    puts "Creating file...\n\n"
    year, month, day, hour, min, sec = Time.now.year, Time.now.month, Time.now.day, Time.now.hour, Time.now.min, Time.now.sec
    day = Time.now.day.to_s.insert(0,"0") if day < 10
    month = Time.now.month.to_s.insert(0,"0") if month < 10
    min = Time.now.min.to_s.insert(0,"0") if min < 10
    sec = Time.now.sec.to_s.insert(0,"0") if sec < 10
    builder = Nokogiri::XML::Builder.new(:encoding => 'UTF-8') do |xml|
      xml.MzIdentML(:id => "",
        :version => "1.1.0",
        'xsi:schemaLocation' => "http://psidev.info/psi/pi/mzIdentML/1.1 ../../schema/mzIdentML1.1.0.xsd",
        'xmlns' => "http://psidev.info/psi/pi/mzIdentML/1.1",
        'xmlns:xsi' => "http://www.w3.org/2001/XMLSchema-instance",
        :creationDate => "#{year}-#{month}-#{day}T#{hour}:#{min}:#{sec}") {
          cvList(xml)
          analysisSoftwareList(xml)
          #analysisSampleCollection(xml)
          provider(xml)  
          auditCollection(xml)
          sequenceCollection(xml)
          analysisCollection(xml)
          analysisProtocolCollection(xml)
          dataCollection(xml)
        }
    end
    
    @pepxml.display_missed_mappings
    File.open(base_file + ".mzid", 'w') {|io| io.puts builder.to_xml}
  end
  
  
  private
  
  # Takes the input filename and the filetype type (tp)
  def base_file
    if @pepxml.type == "pepxml"
      @pepxml.file.chomp('.pep.xml')
    end
  end
  
  # Depth 1. This method name, as well as all the other ones like it, are named after the sections they create.
  def cvList(xml)
    xml.cvList {
      xml.cv(:id => "PSI-MS", :fullName => "Proteomics Standards Initiative Mass Spectrometry Vocabularies", :uri => "http://psidev.cvs.sourceforge.net/viewvc/*checkout*/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo", :version => "2.32.0")
      xml.cv(:id => "UNIMOD", :fullName => "UNIMOD", :uri => "http://www.unimod.org/obo/unimod.obo")
      xml.cv(:id => "UO", :fullName => "UNIT-ONTOLOGY", :uri => "http://obo.cvs.sourceforge.net/*checkout*/obo/obo/ontology/phenotype/unit.obo")
    }
  end
  
  # Depth 1
  def analysisSoftwareList(xml)
  
    #No reference to wich software was used for each SIL, just list the software
    summparams = @pepxml.peptideIdentificationSummaryParams
    
    #combined_software_list = [] 
    #De momento solo tengo X!Tandem. En los 3 SIL. No tiene sentido hacer este combined software list
    #summparams.each do |s| # Software_list Per SIL (AtiO2, Elu1A, Elu2A)
    #  s.search_params_set.searchparams.software_list.each do |soft|
    #   combined_software_list << soft
    #  end
    #end
  
    #combined_software_list.flatten!.uniq! #No software por SIL, sino global!
    softwareList = summparams[0].search_params_set.searchparams.software_list #[0] es el del primer SIL Ati02, me da igual , todos son X!Tandem
    
    xml.AnalysisSoftwareList {
      #combined_software_list.each do |software_name|
      softwareList.each do |software_name, version|  
        xml.AnalysisSoftware(:id => software_name, :version => version) {
          #1 <ContactRole>
          if software_name == "X\\!Tandem"
            xml.ContactRole(:contact_ref => "ORG_GPM" ) {
              xml.Role {
                xml.cvParam(:accession => "MS:1001267", :name => "software vendor",  :cvRef => "PSI-MS")
              }
            }
          elsif software_name == "peptideProphet" or software_name == "proteinProphet"
            xml.ContactRole(:contact_ref => "ORG_ISB" ) {
              xml.Role {
                xml.cvParam(:accession => "MS:1001267", :name => "software vendor",  :cvRef => "PSI-MS")
              }
            }            
          elsif software_name == "pepXML_protXML_2_mzIdentML"
            #ContactRole no es obligatorio : minOccurs:0. Pero creo que es necesario para MIAPE compliance en el validator
            xml.ContactRole(:contact_ref => "ORG_UCM") {
              xml.Role { 
                xml.cvParam(:accession => "MS:1001267", :name => "software vendor",  :cvRef => "PSI-MS")
              }
            }          
          end
          #2 <SoftwareName>
          if @pepxml.findAccession(software_name)[0] != "" 
            xml.SoftwareName {
              array = @pepxml.findAccession(software_name)
              xml.cvParam(:accession => array[0], :name => array[1], :cvRef => "PSI-MS")
            }
          else
            xml.SoftwareName {
              xml.userParam(:name => software_name)
            }
          end
          #3 <Customizations>
          xml.Customizations "k-score scoring algorithm" if software_name == "X\\!Tandem"
        }
      end
    }   
    
  end
  
  
  # Depth 1
  def provider(xml)
    xml.Provider(:analysisSoftware_ref => "pepXML_protXML_2_mzIdentML", :id => "PROVIDER") {
      xml.ContactRole(:contact_ref => "PERSON_DOC_OWNER") {
        xml.Role {
          xml.cvParam(:accession => "MS:1001271", :name => "researcher", :cvRef => "PSI-MS")
        }
      }
    }
  end
  
  #Depth 1
  def auditCollection(xml)
    xml.AuditCollection {
      xml.Person(:id => "PERSON_DOC_OWNER", :firstName => "Vital", :lastName => "Vialas") {
        xml.cvParam(:accession => "MS:1000589", :name => "contact email", :value => "vital@farm.ucm.es", :cvRef => "PSI-MS")
        xml.Affiliation(:organization_ref => "ORG_UCM")
      }
      xml.Organization(:id => "ORG_UCM", :name => "Universidad Complutense de Madrid"){
        xml.cvParam(:accession => "MS:1000588", :name => "contact URL", :value => "http://www.ucm.es", :cvRef => "PSI-MS")
      }
      xml.Organization(:id => "ORG_GPM", :name => "The Global Proteome Machine"){
        xml.cvParam(:accession => "MS:1000588", :name => "contact URL", :value => "http://www.thegpm.org/", :cvRef => "PSI-MS")
      }
      xml.Organization(:id => "ORG_ISB", :name => "Institute for Systems Biology"){
        xml.cvParam(:accession => "MS:1000588", :name => "contact URL", :value => "http://www.systemsbiology.org/", :cvRef => "PSI-MS")
      }
    }
  end
  
  # Depth 1
  def sequenceCollection(xml)
    xml.SequenceCollection {
      dBSequences(xml)
      peptides(xml)
      peptideEvidences(xml)
    }
  end
  
  # Depth 2
  def dBSequences(xml)
    puts "  get proteins from pepxml object, method proteins and write <db_sequence>s "
    proteins = @pepxml.proteins    
    proteins.each do |protein|
      xml.DBSequence(:id => protein[2], :searchDatabase_ref => "SDB_1", :accession => protein[0]) {
        xml.cvParam(:accession => "MS:1001088", :name => "protein description", :cvRef => "PSI-MS", :value => protein[1])
        #Meter Aqui la secuencia!!
        xml.Seq protein[3]
      }
    end
    this_time = Time.now
    puts "done\n"
    puts "Time elapsed: #{this_time - @@beginning_time} seconds\n\n"
  end
  
  # Depth 2
  def peptides(xml)
    puts "get peptides from pepxml object, method peptides and write <peptides>"
    peptides = @pepxml.peptides    
    peptides.each do |peptide|
      pep_id = peptide[0]
      pep_seq, modif_arr = peptide[1][0], peptide[1][1]
      xml.Peptide(:id => pep_id) {
        xml.PeptideSequence  pep_seq 
        if !modif_arr.empty?
          modif_arr.each do |modif|
            cvParam = modif[-1]
            xml.Modification(:residues => modif[0], :location => modif[1], :avgMassDelta => modif[2] ) {
              xml.cvParam(:cvRef => cvParam[0], :accession => cvParam[1], :name => cvParam[2])  if !cvParam.nil?
            }
          end
        end
      }
    end
    this_time = Time.now
    puts "done\n"
    puts "Time elapsed: #{(this_time - @@beginning_time)/60} minutes\n\n"
  end

  # Depth 2
  def peptideEvidences(xml)
    puts "  get results from pepxml object, method results. Write <peptide_evidence>s"
    @pepxml.peptideEvidences.each do |ev|
      xml.PeptideEvidence(
      :dBSequence_ref => "SDB_1_#{ev[4]}",
      #:frame => "",
      :id => "#{ev[0]}_SDB_1_#{ev[4]}_#{ev[5]}_#{ev[6]}", 
      :isDecoy => ev[7],
      :name => "",
      :peptide_ref => ev[0], 
      :start => ev[5],
      :end => ev[6],
      :pre => ev[1],
      :post => ev[2] )
    end

    this_time = Time.now
    puts "done\n"
    puts "Time elapsed: #{(this_time - @@beginning_time) / 60} minutes\n\n"
    
  end

  
  # Depth 1
  def analysisCollection(xml)
    puts "get identification parameters and write 2 <AnalysisCollection>  <AnalysisProtocolCollection>"

    xml.AnalysisCollection {
      #SI_ATiO2 - #SIP_ATiO2 - #SIL_ATiO2 - #SID_ATiO2 
      @summparams = @pepxml.peptideIdentificationSummaryParams
      all_sil_s = []
      @summparams.each do |s|
        sil = s.sil
        spectra_data_ref = sil.sub("SIL", "SID")
        all_sil_s << s.sil
        database = s.search_params_set.searchparams.search_database
        if s.ms_run_params_set.msrunparams.inputraw_path =~ /^\//
          inputraw = s.ms_run_params_set.msrunparams.inputraw_path.split("/")[-1]
        else
          inputraw = s.ms_run_params_set.msrunparams.inputraw_path
        end
        xml.SpectrumIdentification(:id => "SI_#{inputraw}", :spectrumIdentificationProtocol_ref => "SIP_#{inputraw}", :spectrumIdentificationList_ref => sil) {
          xml.InputSpectra(:spectraData_ref => spectra_data_ref)
          xml.SearchDatabaseRef(:searchDatabase_ref => "SDB_1")
        }
      end
      #<ProteinDetection> -maxOccurs = 1 #Static id name is OK
      #<ProteinDetectionList> (under <AnalysisData>) -maxOccurs = 1 
      #<ProteinDetectionProtocol> (under <AnalysisProtocolCollection>) -maxOccurs = 1   TOO ! Static ids! Cool! Easy!
      xml.ProteinDetection(:id => "PD_1", :proteinDetectionList_ref => "PDL_1", :proteinDetectionProtocol_ref => "PDP_1") {
        all_sil_s.each do |sil|
         xml.InputSpectrumIdentifications(:spectrumIdentificationList_ref => "#{sil}")
        end
      }
    }
  end
  
  # Depth 1
  def analysisProtocolCollection(xml)
    xml.AnalysisProtocolCollection {
      SpectrumIdentificationProtocol(xml)
      ProteinDetectionProtocol(xml) #maxOccurs = 1
    }
    this_time = Time.now
    puts "done\n"
    puts "Time elapsed: #{(this_time - @@beginning_time) / 60} minutes\n\n"
  end
  
  # Depth 2
  def SpectrumIdentificationProtocol(xml)
    
    enz_i = 0 # !?! Can't repeat enzyme id ! though I should be able to (trypsin in the 3 different SIPs)
    @summparams.each do |s| #x3 
      if s.ms_run_params_set.msrunparams.inputraw_path =~ /^\//
        inputraw = s.ms_run_params_set.msrunparams.inputraw_path.split("/")[-1]
      else
        inputraw = s.ms_run_params_set.msrunparams.inputraw_path
      end
      searched_mod_param_arrs = s.search_params_set.searchparams.searched_modif_Arr
      additional_search_params = s.search_params_set.searchparams.additional_search_params      
      #search_engine = s.search_params_set.searchparams.software_list[-1][0]
      search_engine = "X\\!Tandem"
      xml.SpectrumIdentificationProtocol(:id => "SIP_#{inputraw}", :analysisSoftware_ref => search_engine) {
        xml.SearchType { #support for ms-ms only #Don't know of any value other than "ms-ms search." Should probably fix this in the future. #MS:1001010 (de novo search), MS:1001031 (spectral library search), MS:1001081 (pmf search), MS:1001584 (combined pmf + ms-ms search)
          xml.cvParam(:accession => "MS:1001083", :name => "ms-ms search", :cvRef => "PSI-MS")
        }     
        xml.AdditionalSearchParams {
          additional_search_params.each do |param|
            xml.userParam(:name => param[0].to_s, :value => param[1].to_s)
          end
        }
        
        xml.ModificationParams {
          searched_mod_param_arrs.each do |m_arr| 
            cvParam = m_arr[4] 
            #puts cvParam          
            xml.SearchModification(:residues => m_arr[0], :massDelta => m_arr[1], :fixedMod => m_arr[3]) {
              xml.cvParam(:cvRef => cvParam[0], :accession => cvParam[1], :name => cvParam[2])
            }
          end    
        }        
        
        xml.Enzymes {
          xml.Enzyme(:id => "Enz_#{enz_i}", :cTermGain => "OH", :nTermGain => "H", :semiSpecific => "0"){
            xml.SiteRegexp  "\<\![CDATA[(\?\<=[KR])(\?\!P)]]\>" 
            xml.EnzymeName{
              xml.cvParam(:accession => "MS:1001251", :name => "Trypsin", :cvRef => "PSI-MS")            
            }

          }        
        }
        enz_i += 1
        
        xml.FragmentTolerance {
          frag_tol = s.search_params_set.searchparams.frag_tol
          xml.cvParam(:accession => "MS:1001412", :name => "search tolerance plus value", :value => frag_tol, :cvRef => "PSI-MS")
          xml.cvParam(:accession => "MS:1001413", :name => "search tolerance minus value", :value => frag_tol, :cvRef => "PSI-MS")
        }
        
        xml.ParentTolerance {
          parent_tol = s.search_params_set.searchparams.parent_tol   
          xml.cvParam(:accession => "MS:1001412", :name => "search tolerance plus value", :value => parent_tol, :cvRef => "PSI-MS")
          xml.cvParam(:accession => "MS:1001413", :name => "search tolerance minus value", :value => parent_tol, :cvRef => "PSI-MS")
        }

        xml.Threshold {
          if @pepxml.threshold == 0
            xml.cvParam(:accession => "MS:1001494", :name => "no threshold", :cvRef => "PSI-MS")
          else
            xml.cvParam(:accession => "MS:?", :name => "?", :cvRef => "PSI-MS", :value => @pepxml.threshold)
          end
        }
            
      }
    end
  end

  # Depth 2
  def ProteinDetectionProtocol(xml)
    #ProteinDetectionProtocol -maxOccurs = 1
    xml.ProteinDetectionProtocol(:id => "PDP_1", :analysisSoftware_ref => "X\\!Tandem" ) {
      xml.Threshold {
        if @protxml.file =~ /highprob/
          xml.userParam(:name => "highprob prot.xml")
        end
      }
    }

  end
  
  # Depth 1
  def dataCollection(xml)
    xml.DataCollection {
      inputs(xml)
      xml.AnalysisData {
        spectrumIdentificationList(xml)
        proteinDetectionList(xml)
      }
    }
  end
  
  # Depth 2
  def inputs(xml)
    xml.Inputs { #The database searched and the input file converted to mzIdentML
      xml.SourceFile( :id => "pep_xml", :location => @pepxml.file )
      xml.SourceFile( :id => "prot_xml", :location => @protxml.file )
      #SearchDatabase x search
      xml.SearchDatabase(:location => @pepxml.database, :id => "SDB_1") {
        xml.DatabaseName {
          xml.userParam(:name => File.basename(@pepxml.database))
        }
      }
      
      @summparams.each do |s|      
        manufacturer = s.ms_run_params_set.msrunparams.manufacturer    
        spectra_data_ref = s.sil.sub("SIL", "SID")
        #<SpectraData> -minOccurs 1 
        #location = s.search_params_set.searchparams.frag_tol
        location = s.ms_run_params_set.msrunparams.inputraw_path + s.ms_run_params_set.msrunparams.raw_data_file_type
        #MALLL!! :location no es eso! tengo que poner la location del mzML !!
        xml.SpectraData(:id => spectra_data_ref, :location => location ) {
          xml.SpectrumIDFormat {
            if manufacturer.include? "Thermo" or manufacturer.include? "thermo" 
              xml.cvParam(:accession => "MS:1000768", :name => "Thermo nativeID format", :cvRef => "PSI-MS")
            end
          }        
        }        
      end
    }
  end


  # Depth 3
  def spectrumIdentificationList(xml)
    @pepxml.peptideIdentificationSummaryParams.each do |s|
      sil = s.sil
      xml.SpectrumIdentificationList(:id => sil) {
        spectrumIdentificationResult(xml, sil)
      }
    end

  end
  
  # Depth 4
  def spectrumIdentificationResult(xml, sil)
  
    puts "writing <spectrum_identification_results> for sil: #{sil}"
    results = @pepxml.results[sil]
    i = 1
    spectra_data_ref = sil.sub("SIL", "SID")
    results.each do |result|
      xml.SpectrumIdentificationResult(:id => "SIR_#{i}_#{sil}", :name => result.name, :spectrumID => "index=#{result.index}", :spectraData_ref => spectra_data_ref) {
        result.items.each do |item|
          ident = item.ident
          #siiID = "SII_#{i}_#{ident.id}"
          siiID = ident.sii_id
          
          #Depth 5
          xml.SpectrumIdentificationItem(
            :id => siiID,
            :calculatedMassToCharge => '%.8f' % ident.mass,  #The 8 dec format is arbitrary. I just felt like it.
            :chargeState => ident.charge,
            :experimentalMassToCharge => '%.8f' % ident.experi,
            :peptide_ref => ident.pep,
            :rank => ident.rank,
            :passThreshold => ident.pass) {
              spectrumIdentificationItemVals(xml, item, siiID)
            }
        end
      }  if result.items.length > 0  #Schema says that SpectrumIdentificationResult can't be empty
      i += 1
    end
    this_time = Time.now
    puts "done\n"
    puts "Time elapsed: #{(this_time - @@beginning_time) / 60} minutes\n\n"
  end
  
  # Depth 6
  def spectrumIdentificationItemVals(xml, item, siiID)
    #OJO! Porque Peptide_ref(s) pueden ser varios por cada SII !!
    pep_id = item.ident.pep
    
    pep_evidences = item.pepEvidences
    
    pep_ev_refs = []
    pep_evidences.each do |ev|
      pep_ev_refs << "#{ev[0]}_SDB_1_#{ev[4]}_#{ev[5]}_#{ev[6]}"
    end
        
    pep_ev_refs.each do |ev|
      xml.PeptideEvidenceRef(:peptideEvidence_ref => ev) 
    end

    cv_params_set, user_params_set = item.vals[0], item.vals[1]
    cv_params_set.each do |cvparam|
      xml.cvParam(:accession => cvparam[0], :name => cvparam[1], :cvRef => "PSI-MS", :value => cvparam[2])
    end
    if !user_params_set.empty?
      user_params_set.each do |userparam|
        xml.userParam(:name => userparam[0], :value => userparam[1]) unless userparam.nil?
      end
    end
   
   #OJORR! En estos cvParams van los scores pero también otras cosas como number of matched / unmatched ions!!
  end
  

  def proteinDetectionList(xml)

    xml.ProteinDetectionList(:id => "PDL_1"){
     puts "getting protein groups from prot_xml object and writing <protein_hypothesis>..."
      protein_groups = @protxml.proteinGroups
      protein_groups.each do |group|
        pag_name = group.name
        xml.ProteinAmbiguityGroup(:id => pag_name)  {
          group.protein_hypothesis_set.each do |protein_hyp|
            proteinHypothesis(xml, protein_hyp)
          end
        }
      end
    }
    this_time = Time.now
    puts "done\n"
    puts "Time elapsed: #{(this_time - @@beginning_time)/60} minutes\n\n"
  
  end  


  def proteinHypothesis(xml, protein_hyp)
   
    #1.upto protein_hyp.n_indist_proteins do
      xml.ProteinDetectionHypothesis(:id => protein_hyp.name, :passThreshold => protein_hyp.pass_threshold) {
        params = protein_hyp.vals

        #El problema es que si tengo indist_proteins > 1 para las siguientes a la primera solo tengo protein_name
      
        protein_hyp.pep_hypothesis_set.each do |pep|
          pep.pep_hypothesis_hash.each do |pep_ev, sii_s|
            xml.PeptideHypothesis(:peptideEvidence_ref => pep_ev) {
              sii_s.each do |sii_id|
                xml.SpectrumIdentificationItemRef(:spectrumIdentificationItem_ref => sii_id)
              end
            }
          end
        end
        #cvParam id: MS:1001093 name: sequence coverage
        if  !params[3] == ""
          xml.cvParam(:accession => "MS:1001093", :name => "sequence coverage", :cvRef => "PSI-MS", :value => params[3]) 
        else 
          xml.cvParam(:accession => "MS:1001093", :name => "sequence coverage", :cvRef => "PSI-MS", :value => "0") 
        end

        xml.userParam(:name => "n_indist_proteins", :value => params[1]) unless params[1] == ""
        xml.userParam(:name => "probability", :value => params[2]) unless params[2] == ""

        xml.userParam(:name => "unique_stripped_peptides", :value => params[4]) unless params[4] == ""
        xml.userParam(:name => "group_sibling_id", :value => params[5]) unless params[5] == ""
        xml.userParam(:name => "total_number_peptides", :value => params[6]) unless params[6] == ""
        xml.userParam(:name => "pct_spectrum_ids", :value => params[7]) unless params[7] == ""
        xml.userParam(:name => "confidence", :value => params[8]) unless params[8] == ""
        xml.userParam(:name => "prot_length", :value => params[9]) unless params[9] == ""
      }        
        
    #end
      
  end
     

end

      #~ protein_hyp.pep_hypothesis_set.each do |pep|
        #~ pep.pep_evs.each do |pep_ev|
          #~ xml.PeptideHypothesis(:peptideEvidence_ref => pep_ev) {
            #~ pep.sii_s.each do |sii_id| 
              #~ xml.SpectrumIdentificationItemRef(:spectrumIdentificationItem_ref => sii_id) 
              #~ #NORR!! <SpectrumIdentificationItemRef> no tiene subelemento userParam !!
              #~ #Quiza todo esto ponerlo mejor en <SpectrumIdentificationItem>
              #~ #params = pep.vals
              #~ #xml.userParam(:name => "initial_probability", :value => params[0]) unless params[0] == ""
              #~ #xml.userParam(:name => "nsp_adjusted_probability", :value => params[1]) unless params[1] == ""
              #~ #xml.userParam(:name => "weight", :value => params[2]) unless params[2] == ""
              #~ #xml.userParam(:name => "is_nondegenerate_evidence", :value => params[3]) unless params[3] == ""
              #~ #xml.userParam(:name => "n_enzymatic_termini", :value => params[4]) unless params[4] == ""
              #~ #xml.userParam(:name => "n_sibling_peptides", :value => params[5]) unless params[5] == ""
              #~ #xml.userParam(:name => "n_sibling_peptides_bin", :value => params[6]) unless params[6] == ""
              #~ #xml.userParam(:name => "n_instances", :value => params[7]) unless params[7] == ""
              #~ #xml.userParam(:name => "exp_tot_instances", :value => params[8]) unless params[8] == ""
              #~ #xml.userParam(:name => "is_contributing_evidence", :value => params[9]) unless params[9] == ""
            #~ end
          #~ }
        #~ end
        #~ 
      #~ end
      
