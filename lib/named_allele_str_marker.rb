require File.expand_path('str_marker.rb', File.dirname(__FILE__))
class NamedAlleleSTRMarker < STRMarker
  ADD_ATTRIBUTES = [:alleles]
  attr_reader *(ADD_ATTRIBUTES)

  def initialize(_offset, _divider, _repeat_formula, _chromosome, _alleles)
    super(_offset, _divider, _repeat_formula, _chromosome)
    @alleles = {'*' => '?'}.update(_alleles)
  end

  def call_allele(sequence)
    alleles[sequence.length.to_s] || alleles['*']
  end

  def self.attributes
    ATTRIBUTES + ADD_ATTRIBUTES
  end
end
