ProteinDetectionParams = Struct.new(:inputraw, :searchEngine, :manufacturer, :model, :ionization, :analyzer, :detector, :raw_data_file_type)

# A simple class that holds a slightly complex structure of data for mzIdentML
#
# @author  ()

#~ class ProteinDetectionParams
#~ 
  #~ attr_accessor :sil, :ms_run_params_set, :search_params_set
#~ 
  #~ def initialize(sil, ms_run_params_set, search_params_set)
    #~ @sil = sil
    #~ @ms_run_params_set = ms_run_params_set
    #~ @search_params_set = search_params_set
  #~ end
#~ 
#~ end


class ProteinDetectionParamsSet

  attr_accessor :msrunparams

  def initialize(ProteinDetectionParams)
    #~ @msrunparams = msrunparams
  end

end


#~ class SearchParamsSet
#~ 
  #~ attr_accessor :searchparams
#~ 
  #~ def initialize(searchparams)
    #~ @searchparams = searchparams  
  #~ end
#~ 
#~ end




