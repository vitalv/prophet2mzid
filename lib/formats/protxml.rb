
require "formats/format.rb"
require "formats/pepxml.rb"
require "protein_group.rb"

# protXML implementation of Format
# @author 2 Vital Vialas  (May 2012)

#class ProtXML < Format
class ProtXML < PepXML
  # Since the prot.xml file is created from the pepxml file and it uses references to SpectrumIdentificationItem and peptides, it makes sense to make this class inherit from PepXML
  # @param [String] file the location of the protXML file
  # @param [String] database the location of the FASTA database that was used by the search engine
	
  #Vital:
  attr_reader :file, :database 
	
  #@protxml_doc = Nokogiri::XML(File.open("SILAC_phos_OrbitrapVelos_1_interact-highprob.prot.xml"))
  # hyph_protxml_obj = ProtXML.new("Hyphal_extract_OrbitrapVelos/XTK_Calb1_0.05_0.05_yes/hyphal_interact-highprob.prot.xml", "C_albicans_SC5314_version_A21-s02-m01-r07_orf_trans_all_plus_contaminants_DECOY.fasta")
  def initialize(file, database)
    #super #ojo aqui al hacer super, estoy invocando initialize en pepxml con lo que se repite el initialize!!
    @file, @database = file, database
    @type = "protxml"
    @doc = Nokogiri::XML(File.open(file))
    @xmlns = ""
    
    # Nokogiri won't parse out the information of an XML file that uses namespaces unless you add xmlns, and vice versa.
    @xmlns = "xmlns:" if hasNamespace
    
  end



  def proteinDetectionParams #Recoge los datos para <AnalysisCollection><ProteinDetection> y para <AnalysisProtocolCollection><ProteinDetectionProtocol>
  
  #<AnalysisCollection> #OJO! En teoría es un ProteinDetection por cada SIL !! Pero creo que en el protXML, estan todas las proteínas detectadas sin referencia a que SIL corresponden
  #<ProteinDetection id="PD_1" proteinDetectionProtocol_ref="PDP_MascotParser_1" proteinDetectionList_ref="PDL_1" activityDate="2011-03-25T13:33:51">
    #<InputSpectrumIdentifications spectrumIdentificationList_ref="SIL_1"/> #minOccurs = 1
  #</ProteinDetection>
  
  #<AnalysisProtocolCollection>
  
   # <ProteinDetectionProtocol id="PDP_MascotParser_1" analysisSoftware_ref="AS_mascot_parser">
     # <AnalysisParams>
      # <cvParam accession="MS:1001316" name="mascot:SigThreshold" cvRef="PSI-MS" value="0.05"/>
      # <cvParam accession="MS:1001317" name="mascot:MaxProteinHits" cvRef="PSI-MS" value="Auto"/>
      # <userParam name="proteinprophet:param" value="0.01"
      # ...
   # </ProteinDetectionProtocol>
   
   # ms_run_params_set = MsRunParamsSet.new(Msrunparams.new(inputraw, searchEngine, manufacturer, model, ionization, analyzer, detector, raw_data_file_type))
    

    prot_summ_header = @doc.xpath("//xmlns:protein_summary_header")
    ref_db = prot_summ_heaader.xpath("./@reference_database").to_s
    source_files= protein_summary_header.xpath("./@source_files").to_s
    min_pep_prob = protein_summary_header.xpath("./@min_peptide_probability").to_s
    min_peptide_weight = protein_summary_header.xpath("./@min_peptide_weight").to_s
    num_predicted_correct_prots = protein_summary_header.xpath("./@num_predicted_correct_prots").to_s
    num_input_1_spectra = protein_summary_header.xpath("./@num_input_1_spectra").to_s
    num_input_2_spectra = protein_summary_header.xpath("./@num_input_2_spectra").to_s
    num_input_3_spectra = protein_summary_header.xpath("./@num_input_3_spectra").to_s
    num_input_4_spectra = protein_summary_header.xpath("./@num_input_4_spectra").to_s
    num_input_5_spectra = protein_summary_header.xpath("./@num_input_5_spectra").to_s
    initial_min_peptide_prob = protein_summary_header.xpath("./@initial_min_peptide_prob").to_s
    total_no_spectrum_ids  = protein_summary_header.xpath("./@total_no_spectrum_ids").to_s
    sample_enzyme = protein_summary_header.xpath("./@sample_enzyme").to_s

    #~ protein_detection_params = 

  
  
    #~ ProteinDetectionParamsSet.new(
  
  end
  
  


  

  def proteinGroups #OJO, igual que para el caso de ProteinDetection y ProteinDetectionProtocol tndria que hacer varias listas, pero no puedo obtener el SIL de referencia!!
  
    all_protein_groups = @doc.xpath("//#{@xmlns}protein_group")
    
    pags = []
    all_protein_groups.each do |group| #uso all_protein_groups[445] porque tiene varias <protein y una de ellas tien 3 indisti_proteins
      group_number = group.xpath("./@group_number").to_s
      group_name = "PAG_#{group_number}"
      
      protein_hypothesis_set = []
      #protein_hypothesis_set = Array de objetos ProteinHypothesis 
      #ProteinHypothesis contiene protein_name, prot_hypothesis_params, threshold, y peptide_hypothesis_set -> Array de objetos PeptideHypothesis
            
      group.xpath(".//#{@xmlns}protein|.//#{@xmlns}indistinguishable_protein").each do |protein|
        protein_name = protein.xpath("./@protein_name").to_s
        protein = protein.parent if protein.node_name == "indistinguishable_protein"
        protein_hypothesis_params = []
        #Estas protein, si no las encuentro en @@prot_acs, no puedo poner el dBSequence_ref en <ProteinHypothesis
        pass_threshold = "true" #eS REQUIRED, "If no such threshold has been set, value of true should be given for all results"
        dbseq_ref = "SDB_1_#{protein_name}"
        n_indist_proteins = protein.xpath("./@n_indistinguishable_proteins").to_s.to_i
        probability = protein.xpath("./@probability").to_s
        percent_coverage = protein.xpath("./@percent_coverage").to_s  #cvParam id: MS:1001093 name: sequence coverage
        unique_stripped_peptides = protein.xpath("./@unique_stripped_peptides").to_s
        group_sibling_id = protein.xpath("./@group_sibling_id").to_s
        total_number_peptides = protein.xpath("./@total_number_peptides").to_s
        pct_spectrum_ids = protein.xpath("./@pct_spectrum_ids").to_s
        confidence = protein.xpath("./@confidence").to_s
        if protein.xpath(".//#{@xmlns}parameter/@name").to_s == "prot_length"
          prot_length = protein.xpath(".//#{@xmlns}parameter/@value").to_s 
        else
          prot_length = ""
        end     
        prot_hypothesis_params = [dbseq_ref, n_indist_proteins, probability, percent_coverage, unique_stripped_peptides, group_sibling_id, total_number_peptides, pct_spectrum_ids, confidence, prot_length]
        
        peptide_hypothesis_set = []        
        peptides = protein.xpath(".//#{@xmlns}peptide") 
        #@@peps (pepxml.rb) [["peptide_223_6", ["CFTAGTNTVTFNDGGK", [["C", 1, "57.0215", ["UNIMOD", "UNIMOD:4", "Carbamidomethyl"]]]]]]       
        
        #OJO! CuIDAO CON ESTO:
        #peptide_207_1 -> orf19.2020, orf19.2021 y orf19.2023
 
        pep_ids = []
        peptides.each do |pep| #Los <peptide > de una <protein>
          #Obtener peptide_id de cada uno de estos <peptide> (peptide_id: peptide_439_6")
          protxml_pepseq = pep.xpath("./@peptide_sequence").to_s
          protxml_pepmod = pep.xpath(".//#{@xmlns}modification_info")
          @@peps.each do |pepxml_pep|
            #SI PARA ESTE protxml_pepseq no encuentro su pepxml_pepseq correspondiente NO entrará en este PAG, porque no podré aportar un <PeptideHypothesis con su peptideEvidence_ref !!!
            pepxml_pepseq = pepxml_pep[1][0]
            pepxml_pepmod = pepxml_pep[1][1]
            if protxml_pepseq == pepxml_pepseq 
              if (protxml_pepmod.empty? and pepxml_pepmod.empty?) or !protxml_pepmod.empty? #En el prot.xml, Si un pep "SPSS[..]NLSK" está sin modificar (peptide_46_7) y modificado (peptide_20_10) aparece bajo un solo elemento <peptide> !
                pep_ids << pepxml_pep[0] 
              end
            end
          end
        end
        
        #Y luego para cada pep_id saco su(s) pep_Evidence(s) y su(s) sii
        pep_ids.each do |pep_id| #estos son los peptidos de cada proteína
          peptide_hypothesis = {} #peptide_hypothesis["peptide_151_10_SDB_1_orf19.5629_45_61"] = ["SIL_Elu2A_SII_589_1"
          @@peptide_evidences.each do |ev| #["peptide_1_1", "F", "K", "0", "orf19.1086", 554, 569]
          
          
          
            if pep_id == ev[0] and protein_name == ev[4]
            
            
            
              pep_ev = "#{ev[0]}_SDB_1_#{ev[4]}_#{ev[5]}_#{ev[6]}"
              sii_s = []
              @@all_items.each do |sii|
                sii_s << sii.ident.sii_id if sii.ident.pep == pep_id
              end
              #peptide_hypothesis_set << PeptideHypothesis.new(pep_evs, sii_s)
              peptide_hypothesis[pep_ev] = sii_s
            end            
          end
          peptide_hypothesis_set << PeptideHypothesis.new(peptide_hypothesis)
        end       
   
#   id: MS:1001098 name: confident distinct peptide sequences cvParam => Lo cojo del numero de PeptideHypothesis
   
        if !peptide_hypothesis_set.empty?
          protein_hypothesis_set << ProteinHypothesis.new(protein_name, n_indist_proteins, pass_threshold, prot_hypothesis_params, peptide_hypothesis_set)
        end
        
      end
  
      pags << ProteinGroup.new(group_name, protein_hypothesis_set) unless protein_hypothesis_set.empty?
      
    end
    
    
    pags
    
          #~ initial_probability = pep.xpath("./@initial_probability").to_s
          #~ nsp_adjusted_probability = pep.xpath("./@nsp_adjusted_probability").to_s
          #~ weight = pep.xpath("./@weight").to_s
          #~ is_nondegenerate_evidence = pep.xpath("./@is_nondegenerate_evidence").to_s
          #~ n_enzymatic_termini = pep.xpath("./@n_enzymatic_termini").to_s
          #~ n_sibling_peptides = pep.xpath("./@n_sibling_peptides").to_s
          #~ n_sibling_peptides_bin = pep.xpath("./@n_sibling_peptides_bin").to_s
          #~ n_instances = pep.xpath("./@n_instances").to_s
          #~ exp_tot_instances = pep.xpath("./@exp_tot_instances").to_s
          #~ is_contributing_evidence = pep.xpath("./@is_contributing_evidence").to_s
          #~ pep_vals = [initial_probability, nsp_adjusted_probability, weight, is_nondegenerate_evidence, is_nondegenerate_evidence, n_enzymatic_termini, n_sibling_peptides, n_sibling_peptides_bin, n_instances, exp_tot_instances, is_contributing_evidence]
        
      
  
  end

  



end





