require File.expand_path('json_able.rb', File.dirname(__FILE__))
class SNPMarker < JSONAble
  ATTRIBUTES = [:position, :snp_type, :reference_sequence, :chromosome]
  attr_reader *(ATTRIBUTES)
#  attr_writer :chromosome

  def initialize(_position, _snp_type, _reference_sequence = nil, _chromosome = 'A')
    @position = _position.to_i
    @snp_type = _snp_type.to_i
    @reference_sequence = _reference_sequence
    @chromosome = _chromosome
  end

  def call_allele(sequence)
    sequence[position]
  end

  def match(sequence)
    [:match, sequence]
  end

  def self.attributes
    ATTRIBUTES
  end

  def is_str?
    false
  end

  def is_snp?
    true
  end
end
