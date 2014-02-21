
Ident = Struct.new(:sii_id, :mass, :charge, :experi, :pep, :rank, :pass)
#~ PepEvidence = Struct.new(:id, :start, :end, :pre, :post, :missedCleavages, :isDecoy, :DBSequence_Ref)

# A simple class that holds a slightly complex structure of data for mzIdentML
#
# @author Jesse Jashinsky (Aug 2010)
class SpectIdResult

  attr_accessor :index, :items, :name

  def initialize(index, name, items)
    @index = index
    @items = items
    @name = name 
  end


end

# A simple class that holds data for SpectIdResult
#
# @author Jesse Jashinsky (Aug 2010)
class SpectIdItem

   
  def initialize(ident)
    @ident = ident
    @pepEvs = nil
    @pepMod = nil
    @vals = []
  end
  
  def ident
    @ident
  end
  
  def pepEvidences
    @pepEvs
  end
 
  def pepEvidences=(pepEvs)
    @pepEvs = pepEvs
  end

  def vals
    @vals
  end
  
  def vals=(vals)
    @vals = vals
  end
  
end
