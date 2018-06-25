require File.expand_path('sequence_helper.rb', File.dirname(__FILE__))
class Primer
  attr_reader :length, :sequence, :first_half, :second_half, :original_sequence, :original_sequence_as_bytes, :first_half_length, :second_half_length, :original_first_half, :original_second_half

  def initialize(_sequence)
    @original_sequence = _sequence
    @original_sequence.freeze
    @original_sequence_as_bytes = _sequence.bytes
    @original_sequence_as_bytes.freeze

    @length = _sequence.length
    @sequence = SequenceHelper::string_or_iupac_regexp(_sequence)
    @sequence.freeze

    @original_first_half = _sequence[0 ... (@length / 2)]
    @first_half_length = @original_first_half.nil? ? 0 : @original_first_half.length
    @first_half = SequenceHelper::string_or_iupac_regexp(@original_first_half)
    @first_half.freeze

    @original_second_half = _sequence[(@length / 2) .. -1]
    @second_half_length = @original_second_half.nil? ? 0 : @original_second_half.length
    @second_half = SequenceHelper::string_or_iupac_regexp(@original_second_half)
    @second_half.freeze
  end

  def empty?
    @length == 0
  end
end
