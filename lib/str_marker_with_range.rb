require File.expand_path('str_marker.rb', File.dirname(__FILE__))
require File.expand_path('sequence_helper.rb', File.dirname(__FILE__))
class STRMarkerWithRange < STRMarker
  RANGE_ATTRIBUTES = [:forward_range, :reverse_range]
  attr_reader *(RANGE_ATTRIBUTES)

  def initialize(_offset, _divider, _repeat_formula, _chromosome, _forward_range = nil, _reverse_range = nil)
    super(_offset, _divider, _repeat_formula, _chromosome)
    @forward_range = _forward_range.to_i
    @reverse_range = _reverse_range.nil? || _reverse_range.to_i == 0 ? -1 : _reverse_range.to_i
  end

  def match(sequence)
    [:match, sequence[@forward_range..@reverse_range]]
  end

  def self.attributes
    ATTRIBUTES + RANGE_ATTRIBUTES
  end
end
