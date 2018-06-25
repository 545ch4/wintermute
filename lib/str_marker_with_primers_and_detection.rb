require File.expand_path('str_marker.rb', File.dirname(__FILE__))
require File.expand_path('sequence_helper.rb', File.dirname(__FILE__))
class STRMarkerWithPrimersAndDetection < STRMarker
  PRIMERS_AND_DETECTION_ATTRIBUTES = [:forward_primer, :reverse_primer, :detection]
  attr_reader *(PRIMERS_AND_DETECTION_ATTRIBUTES)

  def initialize(_offset, _divider, _repeat_formula, _chromosome, forward_primer_sequence = nil, reverse_primer_sequence = nil, _detection = nil)
    super(_offset, _divider, _repeat_formula, _chromosome)
    @forward_primer = forward_primer_sequence
    @reverse_primer = reverse_primer_sequence
    @target = Target.new(forward_primer_sequence.to_s, reverse_primer_sequence.to_s, Target::FORWARD)
    @algo = FastTargetMatchingAlgorithm.new(1.5, nil, {:scan_for_r2 => false})
    @detection = _detection
    @detect_method = lambda do |sequence|
      puts "I don't know what to do!"
      false
    end
    if (_detection['detector'] || _detection[:detector]) && (_detection['values'] || _detection[:values])
      detector = _detection['detector'] || _detection[:detector]
      if detector.to_sym == :does_contain
        @detector_value = Primer.new((_detection['values'] || _detection[:values]).first)
        @detect_method = lambda do |sequence|
          !sequence.index(@detector_value.first_half).nil? || !sequence.index(@detector_value.second_half).nil?
        end
      elsif detector.to_sym == :does_not_contain
        @detector_value = Primer.new((_detection['values'] || _detection[:values]).first)
        @detect_method = lambda do |sequence|
          sequence.index(@detector_value.first_half).nil? && sequence.index(@detector_value.second_half).nil?
        end
      else
        @detect_method = lambda do |sequence|
          puts "I don't know what to do with #{detector.inspect} and #{@detector_values.inspect}!"
          false
        end
      end
    end
  end

  def match(sequence)
    if @detect_method.call(sequence)
      (reason, matching_sequence) = @algo.match(@target, sequence)
      matching_sequence.nil? ? [:str_primers_mismatch, nil] : [:str_primers_match, matching_sequence]
    else
      [:str_detection_fail, nil]
    end
  end

  def self.attributes
    ATTRIBUTES + PRIMERS_AND_DETECTION_ATTRIBUTES
  end
end
