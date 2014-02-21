
#~ Ident = Struct.new(:id, :mass, :charge, :experi, :pep, :rank, :pass)
#~ PepEvidence = Struct.new(:id, :start, :end, :pre, :post, :missedCleavages, :isDecoy, :DBSequence_Ref)

# A simple class that holds a slightly complex structure of data for mzIdentML
#
# @author 

class ProteinGroup

  attr_accessor :name, :protein_hypothesis_set

  def initialize(name, protein_hypothesis_set)
    @name = name
    @protein_hypothesis_set = protein_hypothesis_set
  end
  
end

# A simple class that holds data for 
#
class ProteinHypothesis

  attr_accessor :name, :n_indist_proteins, :pass_threshold, :vals, :pep_hypothesis_set 

  def initialize(name, n_indist_proteins, pass_threshold, vals, pep_hypothesis_set)
    @pass_threshold = pass_threshold
    @name = name
    @n_indist_proteins = n_indist_proteins
    @vals = vals
    @pep_hypothesis_set = pep_hypothesis_set
  end

  #~ def peptideHypothesisSet
    #~ @peptideHypothesisSet
  #~ end
   #~ 
   
  # def peptideHypothesis=(peptideHypothesisSet)
    # @peptideHypothesisSet = peptideHypothesisSet
  # end
  
  
end

class PeptideHypothesis

 # attr_accessor :pep_evs, :sii_s

  #def initialize(pep_evs, sii_s)
   # @pep_evs = pep_evs
    #@sii_s = sii_s
    ##@vals = vals
  #end
  
  attr_accessor :pep_hypothesis_hash
  
  def initialize(pep_hypothesis)
     @pep_hypothesis_hash = pep_hypothesis
  end

end


