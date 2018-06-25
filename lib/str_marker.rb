require File.expand_path('json_able.rb', File.dirname(__FILE__))
class STRMarker < JSONAble
  ATTRIBUTES = [:offset, :divider, :repeat_formula, :chromosome]
  attr_reader *(ATTRIBUTES)
#  attr_writer :chromosome

  def initialize(_offset, _divider, _repeat_formula, _chromosome = 'A')
    @offset = _offset.to_i
    @divider = _divider.to_i
    @repeat_formula = _repeat_formula
    @chromosome = _chromosome.to_s
  end

  def call_allele(sequence)
    ('%i.%i'%([(sequence.length - offset) / divider, (sequence.length - offset) % divider])).sub('.0', '')
  end

  def match(sequence)
    [:match, sequence]
  end

  def self.attributes
    ATTRIBUTES
  end

  def is_str?
    true
  end

  def is_snp?
    false
  end
end
