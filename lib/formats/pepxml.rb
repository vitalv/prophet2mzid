#require "formats/format.rb"
require "#{File.expand_path(File.dirname(__FILE__))}/format.rb" 
#require "natcmp.rb"
require "#{File.expand_path(File.dirname(__FILE__))}/../natcmp.rb" 
require "ms/fasta.rb"
#require "peptide_identification_summary_params.rb"
require "#{File.expand_path(File.dirname(__FILE__))}/../peptide_identification_summary_params.rb" 

# The pepXML implementation of Format
#
# @author Jesse Jashinsky (Aug 2010)
# @author 2 Vital Vialas  (May 2012)
class PepXML < Format
  # @param [String] file the location of the pepXML file
  # @param [String] database the location of the FASTA database that was used by the search engine

  #Vital:
  attr_reader :file, :database, :databasename, :type, :date, :threshold, :numsequences, :searchEngine

   #my_pepxml_obj = PepXML.new("SILAC_phos_OrbitrapVelos_1_interact-ipro-filtered.pep.xml", "C_albicans_SC5314_version_A21-s02-m01-r07_orf_trans_all_plus_contaminants_DECOY.fasta")
   #hyph_pepxml_obj = PepXML.new("Hyphal_extract_OrbitrapVelos/XTK_Calb1_0.05_0.05_yes/hyphal_interact-ipro-filtered.pep.xml", "C_albicans_SC5314_version_A21-s02-m01-r07_orf_trans_all_plus_contaminants_DECOY.fasta")

   @@all_items = []

  def initialize(file, database)
    super
    @file, @database = file, database
    @type = "pepxml"
    @doc = Nokogiri::XML(File.open(file))
    #@pepxml_doc = Nokogiri::XML(File.open("SILAC_phos_OrbitrapVelos_1_interact-ipro-filtered.pep.xml"))
    @xmlns = ""
    @numsequences = 0
    @proteinIndices = []

    @threshold = 0

    # Nokogiri won't parse out the information of an XML file that uses namespaces unless you add xmlns, and vice versa.
    @xmlns = "xmlns:" if hasNamespace

    findAllPepLocations # Finds all pep locations and puts them in array: [[peptide, protein, start, end]]

    peptideIdentificationSummaryParams
    @searchEngine = ""

    temp = database.split("/")
    @databaseName = temp[temp.length-1]

    puts "pepxml object initialized"

  end


  def proteins

    puts "getting protein names, descriptions and sequences"
    #OJO, ESTAS PROTEINAS, que son las que uso para SequenceCollection en el .mzid, ¿debería cogerlas del prot.xml  ?

    #@pros[1] es  ["orf19.2167", "orf19.2167 CGDID:CAL0005395 Description.......ization", "DBSeq_1_orf19.2167"]
    allHits = @doc.xpath("//#{@xmlns}search_hit/@protein|//#{@xmlns}search_hit/@protein_descr")
    altProts = @doc.xpath("//#{@xmlns}alternative_protein/@protein|//#{@xmlns}alternative_protein/@protein_descr")
    allHits = allHits + altProts

    pros = []
    i = 0

    #Include pep.xml <alternative_protein    para los <ProteinAmbiguityGroup

    while i < allHits.length
      pro = proteinID(allHits[i].to_s)
      if allHits[i+1].name == "protein_descr"
        #pros << [pro, allHits[i+1].to_s, "SDB_1_#{pro}"]
        pros << [pro, allHits[i+1].to_s, "SDB_1_#{pro}"]
        i += 2
      elsif allHits[i+1].name == "protein"
        pros << [pro, "", "SDB_1_#{pro}"]
        i += 1
      end
    end

    @pros = pros.uniq

    @pros.collect { |p| p.push proteinSeq(p[0]) }

    #@@prot_acs = @pros.collect { |p| p[0] }
    #@@dbSeq_refs = @pros.collect {|p| p[2] }

    @pros


  end


  def peptides

    #These peptides will be the unique <Peptide > entries under <SequenceCollection> in the mzid file
    peps = []
    @doc.xpath("//#{@xmlns}search_hit").each do |search_hit|
      modifications_arr = []
      pepseq = search_hit.xpath("./@peptide").to_s
      modifications_arr = getModifications(pepseq, search_hit)
      peps << [pepseq, modifications_arr]
    end

    peps.uniq! 
    #Since these peptides will be the list of every unique sequence under <SequenceCollection> in the .mzid file, get uniq! array

    i = 0
    first = 1
    second = 1
    while i < peps.length
      peps[i] = ["peptide_#{first}_#{second}", peps[i]]
      i += 1
      second += 1
      if second == 11
        first += 1
        second = 1
      end
    end

    @@peps = peps
    peps  # ["peptide_2_6", ["LPAAIDSK", [["L", 1, "6.0201", ["MS", "MS:1001460", "unknown modification"]]]]]

  end



  def peptideEvidences 
    #peptideEvidence is MzIdentMl nomenclature, it is a peptide sequence in a particular protein position
    #this information is under search_hit in protxml
    pep_ev_things = []
    @doc.xpath("//#{@xmlns}search_hit").each do |search_hit|
      pepseq = search_hit.xpath("./@peptide").to_s
      modifications_arr = getModifications(pepseq, search_hit)
      pre = search_hit.xpath("./@peptide_prev_aa").to_s
      post = search_hit.xpath("./@peptide_next_aa").to_s
      missedCleavages = search_hit.xpath("./@num_missed_cleavages").to_s
      pro = proteinID(search_hit.xpath("./@protein").to_s)
      startVal, endVal = pepLocation(search_hit, pro, pepseq)
      pep_ev_things << [pepseq, modifications_arr, pre, post, missedCleavages, pro, startVal, endVal]
      if !search_hit.xpath(".//#{@xmlns}alternative_protein").empty?
        alt_prot_arr = search_hit.xpath(".//#{@xmlns}alternative_protein")
        alt_prot_arr.each do |alt_prot_node|
          pro = proteinID(alt_prot_node.xpath("./@protein").to_s)
          startVal, endVal = pepLocation(alt_prot_node, pro, pepseq)
          pre, post = alt_prot_node.xpath("./@peptide_prev_aa").to_s, alt_prot_node.xpath("./@peptide_next_aa").to_s
          pep_ev_things << [pepseq, modifications_arr, pre, post, missedCleavages, pro, startVal, endVal]
        end
      end
    end
    pep_ev_things.uniq!
    peptide_evidences = []
    @@peps.each do |pep|
      pep_ev_things.each do |pep_ev|
        if pep[1][0] == pep_ev[0] and pep[1][1] == pep_ev[1] #Same peptide (seq+mod). Grab rest of pep_ev_things and check if they are also the same
          peptide_evidences << [pep[0], pep_ev[2], pep_ev[3], pep_ev[4], pep_ev[5], pep_ev[6], pep_ev[7] ]
                             #["pep_1_1",   "K",   "D",        "0",    "orf19.6313.2",   123,     141]
        end
      end
    end
    peptide_evidences.uniq!
    @@peptide_evidences = peptide_evidences
    peptide_evidences

  end


  def peptideIdentificationSummaryParams

    #summary_params son parametros que van en <SpectrumIdentification bajo <AnalysisCollection en el mzid

    summary_params = []
    spect_id_lists = []


    #SOFTWARE LIST MAY BE DIFFERENT IN EACH SIL, BUT I'M INCLUDING SAME LIST IN ALL SILs
    software_list = {} #software_list[software_name] = version
    software_list["pepXML_protXML_2_mzIdentML"] = "1.0"
    software_list["proteinProphet"] = "Insilicos_LabKey_C++ (TPP v4.5 RAPTURE rev 0, Build 201109301446 (linux))"
    software_list["peptideProphet"] = "Insilicos_LabKey_C++ (TPP v4.5 RAPTURE rev 0, Build 201109301446 (linux))"
    software_list["X\\!Tandem"] = "x! tandem 2010.10.01.1 (LabKey, Insilicos and ISB)"
    #software_list = @doc.xpath("//#{@xmlns}analysis_summary/@analysis").select {|a| a.to_s != "database_refresh" }.collect { |a| a.to_s}


    @doc.xpath("//#{@xmlns}msms_run_summary").each do |msmsrun| #This may differ from one msmsrun to another
      inputraw_path = msmsrun.xpath("./@base_name").to_s # x3 (AtiO2, Elu1A, Elu2A)
      inputraw = inputraw_path.split("/")[-1] if inputraw_path =~ /^\//
      sil_id = "SIL_#{inputraw}"
      spect_id_lists << sil_id
      manufacturer = msmsrun.xpath("./@msManufacturer").to_s
      model = msmsrun.xpath("./@msModel").to_s
      ionization = msmsrun.xpath("./@msIonization").to_s
      analyzer = msmsrun.xpath("./@msMassAnalyzer").to_s
      detector = msmsrun.xpath("./@msDetector").to_s
      raw_data_file_type = msmsrun.xpath("./@raw_data").to_s

      #Todo esto tengo que ponerlo en algun sitio en el mzid!! pero donde??
      #Es que ese tipo de parametros (detector, analizador, etc...) es mas de mzml

      ms_run_params_set = MsRunParamsSet.new(Msrunparams.new(inputraw_path, searchEngine, manufacturer, model, ionization, analyzer, detector, raw_data_file_type))

      #enzyme
      #enzyme_specificity

      search_summ = msmsrun.xpath(".//#{@xmlns}search_summary")
      search_database = search_summ.xpath(".//#{@xmlns}search_database/@local_path").to_s
      precursor_mass_type = search_summ.xpath("./@precursor_mass_type").to_s
      fragment_mass_type = search_summ.xpath("./@fragment_mass_type").to_s

      modifications = search_summ.xpath(".//#{@xmlns}aminoacid_modification")
      searched_modif_Arr = []
      if !modifications.empty?
        modifications.each do |mod|
          residue = mod.xpath("./@aminoacid").to_s
          mass_delta = mod.xpath("./@massdiff").to_s
          mass = mod.xpath("./@mass").to_s.to_f
          variable = mod.xpath("./@variable").to_s
          cv_param_arr = getModificationCvParam(residue, mass_delta)
          searched_modif_Arr << [residue, mass_delta, mass, variable2fixed(variable), cv_param_arr]
        end
      end

      parameters = search_summ.xpath(".//#{@xmlns}parameter")
      additionalSearchParamsArr = []
      frag_tol, parent_tol = "", ""
      spectra_file_path = ""
      if !parameters.empty?
        parameters.each do |param|
          name = param.xpath("./@name").to_s
          value = param.xpath("./@value").to_s
          frag_tol = value if name == "spectrum, fragment monoisotopic mass error"
          parent_tol = value if name == "spectrum, parent monoisotopic mass error minus"
          spectra_file_path = value if name == "spectrum, path"
          additionalSearchParamsArr << [name, value]
        end
      end

      #SOFTARE LIST SAME FOR ALL SILs
      #more_software = []
      #search_engine = msmsrun.xpath("./@search_engine").to_s #DONE:  software_list << "X\!Tandem"
      #analysis_software = msmsrun.xpath(".//#{@xmlns}analysis_timestamp/@analysis").select {|a| a.to_s != "database_refresh"}.collect {|a| a.to_s} #peptideprophet, database_refresh; (ONE per <msms_run_summary> )
      #more_software = [search_engine] + analysis_software
      #software_list << more_software


      search_params_set = SearchParamsSet.new(Searchparams.new(search_database, software_list, precursor_mass_type, fragment_mass_type, frag_tol, parent_tol, spectra_file_path, searched_modif_Arr, additionalSearchParamsArr))

      summary_params << SummaryParams.new(sil_id, ms_run_params_set, search_params_set)

    end

    @spect_id_lists = spect_id_lists

    summary_params

  end



  def results

    puts "getting results. search_hits / items "

    result_set = {}
    all_items = []

    @doc.xpath("//#{@xmlns}msms_run_summary").each_with_index do |msmsrun, i| #(x3)
      sil = @spect_id_lists[i]
      queries = msmsrun.xpath(".//#{@xmlns}spectrum_query")
      indicies = msmsrun.xpath(".//#{@xmlns}spectrum_query/@spectrum").collect {|index| index.to_s}
      indicies = indicies.sort {|x,y| String.natcmp(x, y)}
      results = []

      i = 1
      queries.each do |query|
        charge = query.xpath("./@assumed_charge").to_s.to_i
        hits = query.xpath(".//#{@xmlns}search_hit")
        items = []
        rank = 1

        hits.each do |hit|
          items << getItem(sil, hit, i, rank, charge)
          rank += 1
        end
        i += 1

        name = query.xpath("./@spectrum").to_s
        results << SpectIdResult.new(indicies.index(name), name, items)

        all_items << items.flatten

      end

      result_set[sil] = results

    end

    @@all_items = all_items.flatten

    return result_set
    #return @@all_items

  end


  private

  # Checks if the pepXML file used namespaces
  # @return [Boolean] true if it uses namespaces, false if not
  def hasNamespace
    if @doc.xpath("msms_pipeline_analysis").to_s.length == 0
      true
    else
      false
    end
  end

  # @param [Nokogiri] hit the spectra hit information
  # @param [Integer] rank the rank
  # @param [Integer] charge the charge
  # @return [SpectIdItem] the result item
  def getItem(sil, hit, i, rank, charge)
    mass = hit.xpath("./@calc_neutral_pep_mass").to_s.to_f
    diff = hit.xpath("./@massdiff").to_s.to_f
    scores = hit.xpath(".//#{@xmlns}search_score")
    pepseq = hit.xpath("./@peptide").to_s
    pro = hit.xpath("./@protein").to_s
    mod = hit.xpath(".//#{@xmlns}modification_info")
    #mod_count = mod.xpath(".//#{@xmlns}mod_aminoacid_mass").count
    mods = getModifications(pepseq, hit)

    #OJO CON ESTO: (Mirar Spectrum Identification Protocol donde se define cuál es el threshold, by default set to "no threhsold" )
    #el subelemento <PassThreshold> es required en en <SpectrumIdentificationItem>
    #Por ahora estoy poniendo siempre pass = true pero habra que verlo, no??
    #pass_threshold = "true"

    #OJO con este tema: 1 seq, 2 pep (con y sin modif)
    #peptide_46_7 => ["peptide_46_7", ["SPSSDEVVEYGDLNSANNSANLSK", [] ]]
    #peptide_20_10 => ["peptide_20_10", ["SPSSDEVVEYGDLNSANNSANLSK", [ ["L", 13, "6.0201", ["MS", "MS:1001460", "unknown"]] ] ]]

    pep_id = ""
    pepevs_per_pepid = []
    @@peps.each do |thisPep|
      # 1st check peptide_id (seq+mod)
      if pepseq == thisPep[1][0]
        if (mods == thisPep[1][1]) or (mods.empty? and thisPep[1][1].empty? )
          pep_id = thisPep[0]
          # Then get evidences for that peptide_id
          @@peptide_evidences.each do |ev|
            pepevs_per_pepid << ev if ev[0] == pep_id #and pro == ev[4]
          end
        end

      end
    end

    #SpectIdItem is a class in spect_id_result.rb -> (:id, :mass, :charge, :experi, :pep, :rank, :pass)
    sii_id = "#{sil}_SII_#{i}_#{rank}"
    item = SpectIdItem.new(Ident.new(sii_id, calMass(mass, charge), charge, experiMass(mass, charge, diff), pep_id, rank, "true"))

    obo_mapping_scoreArr, non_obo_mapping_scoreArr = [] , []

    searchEngine = hit.parent.parent.parent.xpath("./@search_engine").to_s

    scores.each do |score|
      id, name = findAccession(conformScoreName(score.xpath("./@name").to_s, searchEngine))
      if id != ""
        obo_mapping_scoreArr << [id, name, score.xpath("./@value").to_s]
      else
        non_obo_mapping_scoreArr << [name, score.xpath("./@value").to_s]
      end
    end


    #METER TAMBIéN AQUI peptideprophet probability Y interprophet probability
    prophet_probs = []
    if hit.xpath(".//#{@xmlns}analysis_result")
      hit.xpath(".//#{@xmlns}analysis_result").each do |a|
        if a.xpath("./@analysis").to_s == "peptideprophet"
          name = "peptide_prophet_probability"
          val = a.xpath(".//#{@xmlns}peptideprophet_result/@probability").to_s
          prophet_probs << [name, val]
        elsif a.xpath("./@analysis").to_s == "interprophet"
          name = "interprophet_probability"
          val = a.xpath(".//#{@xmlns}interprophet_result/@probability").to_s
          prophet_probs << [name, val]
        end
      end
    end

    non_obo_mapping_scoreArr = non_obo_mapping_scoreArr + prophet_probs

    item.vals = [obo_mapping_scoreArr] + [non_obo_mapping_scoreArr]


    # <PeptideEvidenceRef> es un subelemento de <SpectrumIdentificationItem> obligatorio
    #En algunos SII puede haber varios <PeptideEvidenceRef>  -> misma secuencia de peptido en dif. proteinas o posiciones dentro de una prot
    item.pepEvidences = pepevs_per_pepid

    item

  end




  def getModifications(pepseq, hit)
    pep_modifications_Arr = []
    modifications = hit.xpath(".//#{@xmlns}modification_info")
    if !modifications.empty?
      residue, location, mass  = "", 0, 0
      modifications.xpath(".//#{@xmlns}mod_aminoacid_mass").each do |mod|
        location = mod.xpath("./@position").to_s.to_i
        residue = pepseq.split(%r{\s*})[location-1]
        mass = mod.xpath("./@mass").to_s
        massdiff = ""
        hit.parent.parent.parent.xpath(".//#{@xmlns}aminoacid_modification").each do |m|
          aa = m.xpath("./@aminoacid").to_s
          mod_mass = m.xpath("./@mass").to_s.to_f.truncate.to_s
          if  aa == residue and mod_mass == mass.to_f.truncate.to_s
            massdiff = m.xpath("./@massdiff").to_s
          end
          #not sure why, but L , in the pep.xml. appears like so: FGDDDS[167]DDEFDHDEL[119] <mod_aminoacid_mass position="15" mass="0" />
          massdiff = "6.0201" if  residue == "L" and mass == "0"
        end
        cvParam = getModificationCvParam(residue, massdiff)
        pep_modifications_Arr << [residue, location, massdiff, cvParam]
      end
    end
    pep_modifications_Arr
  end


  def getModificationCvParam(residue, massdiff)
    #'%.2f' %  "57.0215".to_f
    cvParam =
      case [residue, massdiff]
        when ["C", "57.0215"]
          ["UNIMOD", "UNIMOD:4", "Carbamidomethyl"]
        when ["M", "15.9949"]
          ["UNIMOD", "UNIMOD:35", "Oxidation"]
        when ["L", "6.0201"]
          ["UNIMOD", "UNIMOD:188", "13C(6) Silac label"]
        when ["R", "6.0201"]
          ["UNIMOD", "UNIMOD:188", "13C(6) Silac label"]
        when ["S", "79.9663"]
          ["UNIMOD", "UNIMOD:21", "Phospho"]
        when ["T", "79.9663"]
          ["UNIMOD", "UNIMOD:21", "Phospho"]
        when ["Y", "79.9663"]
          ["UNIMOD", "UNIMOD:21", "Phospho"]
        when ["Q", "-17.0265"]
          ["UNIMOD", "UNIMOD:28", "Gln->pyro-Glu"]
        when ["C", "-17.0265"]
           ["UNIMOD", "UNIMOD:385", "Ammonia-loss"]
        when ["E", "-18.0106"]
           ["UNIMOD", "UNIMOD:27", "Gln->pyro-Glu"]
      end
    cvParam
  end


  # @param [Nokogiri] hit the spectra hit information
  # @param [String] pro the protein
  # @param [String] pep the peptide
  # @return [Integer, Integer] the start and end location of the peptide
  def pepLocation(hit, pro, pep)
    @locations.each do |location|
      if location[0] == pep && location[1] == pro
        return location[2], location[3]
      end
    end
    return 0, 0    #In case it doesn't find anything
  end

  # Finds all peptide locations and puts them in an array in the format: [[peptide, protein, start, end]]
  def findAllPepLocations
    hits = @doc.xpath("//#{@xmlns}search_hit")
    #alt_prots = @doc.xpath("//#{@xmlns}alternative_protein")
    all = []
    @locations = []
    i = 0

    # Parses out each peptide and protein
    hits.each do |hit|
      all << [hit.xpath("./@peptide").to_s, proteinID(hit.xpath("./@protein").to_s)]
      if !hit.xpath(".//#{@xmlns}alternative_protein").empty?
        alt_prot_arr = hit.xpath(".//#{@xmlns}alternative_protein")
        alt_prot_arr.each do |alt_prot_node|
          pro = proteinID(alt_prot_node.xpath("./@protein").to_s)
          all << [hit.xpath("./@peptide").to_s, proteinID(alt_prot_node.xpath("./@protein").to_s)]
        end
      end
      i += 1
    end


    all.uniq!
    dataHash = Hash.new

    @database = @databaseName
    Ms::Fasta.foreach(@database) do |entry|
#     Ms::Fasta.foreach("./spec/test_files/C_albicans_SC5314_version_A21-s02-m01-r07_orf_trans_all_plus_contaminants_DECOY.fasta") do |entry|
      @numsequences += 1
      pID = ""
      if @database == "C_albicans_SC5314_version_A21-s02-m01-r07_orf_trans_all_plus_contaminants_DECOY.fasta" #if de Vital
        candi_orf_header = entry.header.split("\s")[0] #if entry.header =~ /^orf/
        pID = proteinID(candi_orf_header)
      else #esto es como estaba. Quitar el if, este else y el end y dejar solo la linea sig
        pID = proteinID(entry.header)
      end
      dataHash[pID] = entry.sequence
      @proteinIndices << pID
    end

    all.each do |set|
      if dataHash[set[1]] != nil
        startVal = dataHash[set[1]].scan_i(set[0])[0]
        @locations << [set[0], set[1], startVal + 1, startVal + set[0].length] unless startVal == nil
      end
    end
  end



  def proteinSeq(protein)
    @database = @databaseName
    seq = ""
    Ms::Fasta.foreach(@database) do |entry|
      candi_orf_header = entry.header.split("\s")[0] #if entry.header =~ /^orf/
      pID = proteinID(candi_orf_header)
      seq = entry.sequence if protein == pID
    end
    seq
  end



  def proteinID(protein)
    #If a protein ID contains a "|", then it contains more than just the ID
    if protein.include?('|')
      arr = protein.split("|")[1].split(":")
      if arr.length == 1
        arr[0]
      else
        arr[1]
      end
    #If there's no characters, then it's an index. I don't fully understand regexp, but this works.
    elsif (protein =~ /[A-Z]/i) == nil
      @proteinIndices[protein.to_i]
    else
      protein
    end
  end


  def variable2fixed(boolthing)
    fixed =
      case boolthing
        when "N"
          "true"
        when "Y"
          "false"
      end
    fixed
  end


end




# For quickly getting the start and end indexes of a string
class String
  def scan_i seq
    pos = 0
    ndx = []
    slen = seq.length

    while i = index(seq,pos)
      ndx << i
      pos = i + slen
    end

    ndx
  end
end
