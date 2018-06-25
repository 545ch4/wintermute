require File.expand_path('str_marker.rb', File.dirname(__FILE__))
require File.expand_path('sequence_helper.rb', File.dirname(__FILE__))
class STRMarkerWithPrimers < STRMarker
  PRIMERS_ATTRIBUTES = [:forward_primer, :reverse_primer]
  attr_reader *(PRIMERS_ATTRIBUTES)

  def initialize(_offset, _divider, _repeat_formula, _chromosome, forward_primer_sequence = nil, reverse_primer_sequence = nil)
    super(_offset, _divider, _repeat_formula, _chromosome)
    @forward_primer = forward_primer_sequence
    @reverse_primer = reverse_primer_sequence
    @target = Target.new(forward_primer_sequence.to_s, reverse_primer_sequence.to_s, Target::FORWARD)
    @algo = FastTargetMatchingAlgorithm.new(1.5, nil, {:scan_for_r2 => false})
  end

  def match(sequence)
    (reason, matching_sequence) = @algo.match(@target, sequence)
    matching_sequence.nil? ? [:str_primers_mismatch, nil] : [:str_primers_match, matching_sequence]
  end

  def self.attributes
    ATTRIBUTES + PRIMERS_ATTRIBUTES
  end
end
