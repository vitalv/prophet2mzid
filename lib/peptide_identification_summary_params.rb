Msrunparams = Struct.new(:inputraw_path, :searchEngine, :manufacturer, :model, :ionization, :analyzer, :detector, :raw_data_file_type)
Searchparams = Struct.new(:search_database, :software_list, :precursor_mass_type, :fragment_mass_type, :frag_tol, :parent_tol, :spectra_file_path, :searched_modif_Arr, :additional_search_params)

# A simple class that holds a slightly complex structure of data for mzIdentML
#
# @author  ()

class SummaryParams

  attr_accessor :sil, :ms_run_params_set, :search_params_set

  def initialize(sil, ms_run_params_set, search_params_set)
    @sil = sil
    @ms_run_params_set = ms_run_params_set
    @search_params_set = search_params_set
  end

end


class MsRunParamsSet

  attr_accessor :msrunparams

  def initialize(msrunparams)
    @msrunparams = msrunparams
  end

end


class SearchParamsSet

  attr_accessor :searchparams

  def initialize(searchparams)
    @searchparams = searchparams  
  end

end




