require File.expand_path('str_marker.rb', File.dirname(__FILE__))
class STRMarkerWithDetection < STRMarker
  DETECT_ATTRIBUTES = [:detection]
  attr_reader *(DETECT_ATTRIBUTES)

  def initialize(_offset, _divider, _repeat_formula, _chromosome, _detection)
    super(_offset, _divider, _repeat_formula, _chromosome)
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
    @detect_method.call(sequence) ? [:str_detection_success, sequence] : [:str_detection_fail, nil]
  end

  def self.attributes
    ATTRIBUTES + DETECT_ATTRIBUTES
  end
end
