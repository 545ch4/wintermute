require 'highline'

class PairedEnds
  attr_reader :forward, :forward_primer, :reverse, :reverse_primer

  def initialize(_forward_primer, _reverse_primer)
    @forward_primer = _forward_primer
    @forward_primer.freeze

    @reverse_primer = _reverse_primer
    @reverse_primer.freeze
  end

  def forward
    @forward_primer
  end

  def reverse
    @reverse_primer
  end
end

class Alignment
  attr_reader :target, :sequence, :alignment, :score
  def initialize(_target, _alignment, _sequence, _score)
    @target = _target.to_s
    @target.freeze

    @alignment = _alignment.to_s
    @alignment.freeze

    @sequence = _sequence.to_s
    @sequence.freeze

    @score = _score.to_f
    @score.freeze
  end

  def to_s
    "#{@target}\n#{@alignment}\n#{@sequence}\nscore = #{@score}"
  end

  def print
    (columns, lines) = HighLine::SystemExtensions.terminal_size
    counter_max = [@target.length, @alignment.length, @sequence.length].max
    counter_width = (counter_max / columns).ceil.to_s.length
    columns = columns - (counter_width + 1)
    i = 0
    j = 0
    while(i < counter_max)
      puts '%0*i '%([counter_width, j]) + @target[i...(i + columns)].to_s
      puts '%0*i '%([counter_width, j]) + @alignment[i...(i + columns)].to_s
      puts '%0*i '%([counter_width, j]) + @sequence[i...(i + columns)].to_s
      i += columns
      j += 1
    end
    puts "score = #{@score}\n\n"
  end
end

class PairedEndsAlignment < Alignment
  attr_reader :flank_or_repeat

  def initialize(_target, _alignment, _sequence, _score)
    super(_target, _alignment, _sequence, _score)
    @flank_or_repeat = _target[_alignment.index('>').._alignment.rindex('>')].to_s
    @flank_or_repeat.freeze
  end
end

class AlignmentAlgorithm
  MATCHES = {
    'A'.ord => ['A'.ord, 'R'.ord, 'W'.ord, 'M'.ord, 'D'.ord, 'H'.ord, 'V'.ord, 'N'.ord].freeze,
    'T'.ord => ['T'.ord, 'Y'.ord, 'W'.ord, 'K'.ord, 'B'.ord, 'D'.ord, 'H'.ord, 'N'.ord].freeze,
    'G'.ord => ['G'.ord, 'R'.ord, 'S'.ord, 'K'.ord, 'B'.ord, 'D'.ord, 'V'.ord, 'N'.ord].freeze,
    'C'.ord => ['C'.ord, 'Y'.ord, 'S'.ord, 'M'.ord, 'B'.ord, 'H'.ord, 'V'.ord, 'N'.ord].freeze,
    'Y'.ord => ['C'.ord, 'T'.ord, 'Y'.ord, 'B'.ord, 'H'.ord, 'N'.ord].freeze,
    'R'.ord => ['A'.ord, 'G'.ord, 'R'.ord, 'D'.ord, 'V'.ord, 'N'.ord].freeze,
    'S'.ord => ['G'.ord, 'C'.ord, 'S'.ord, 'B'.ord, 'V'.ord, 'N'.ord].freeze,
    'W'.ord => ['A'.ord, 'T'.ord, 'W'.ord, 'D'.ord, 'H'.ord, 'N'.ord].freeze,
    'K'.ord => ['T'.ord, 'G'.ord, 'K'.ord, 'B'.ord, 'D'.ord, 'N'.ord].freeze,
    'M'.ord => ['A'.ord, 'C'.ord, 'M'.ord, 'H'.ord, 'V'.ord, 'N'.ord].freeze,
    'B'.ord => ['C'.ord, 'G'.ord, 'T'.ord, 'B'.ord, 'N'.ord].freeze,
    'D'.ord => ['A'.ord, 'G'.ord, 'T'.ord, 'D'.ord, 'N'.ord].freeze,
    'H'.ord => ['A'.ord, 'C'.ord, 'T'.ord, 'H'.ord, 'N'.ord].freeze,
    'V'.ord => ['A'.ord, 'C'.ord, 'G'.ord, 'V'.ord, 'N'.ord].freeze,
    'N'.ord => ['A'.ord, 'C'.ord, 'G'.ord, 'T'.ord, 'N'.ord, 'Y'.ord, 'R'.ord, 'S'.ord, 'W'.ord, 'K'.ord, 'M'.ord, 'B'.ord, 'D'.ord, 'H'.ord, 'V'.ord].freeze
  }.freeze

  SYMBOLS = {
    match: '|'.freeze,
    gap: '-'.freeze,
    mismatch: '*'.freeze,
    target_gap: 'T'.freeze,
    sequence_gap: 'S'.freeze
  }.freeze


  PENALTIES = {
    match: 0.0,
    mismatch: 0.8,
    gap: 0.9,
    consecutive_gap: 0.4,
    trailing_or_heading_sequence_gap: 0.9,
    consecutive_trailing_or_heading_sequence_gap: 0.0
  }.freeze

  attr_reader :penalties, :max_score, :options

  def initialize(_max_score = 1.0, _penalties = nil)
    @penalties = (_penalties ? _penalties : PENALTIES)
    @max_score = _max_score
    @options = {}
  end

  def raw_align(_target, _sequence)
    target = _target.bytes
    sequence = _sequence.bytes
    max_trailing_target_gaps = 1 + ((@max_score - @penalties[:gap]) / @penalties[:consecutive_gap]).floor
    best_score = @max_score + 0.1
    best_alignment_obj = [nil, best_score]
    ((-max_trailing_target_gaps) .. (target.size - sequence.size + max_trailing_target_gaps)).each do |i|
      if i < 0
        alignment_obj = _align([SYMBOLS[:target_gap] * i.abs, @penalties[:gap] + (i.abs.pred * @penalties[:consecutive_gap])], target, 0, sequence, i.abs, best_score)
      elsif i == 0
        alignment_obj = _align(['', 0.0], target, i, sequence, 0, best_score)
      else
        alignment_obj = _align([SYMBOLS[:sequence_gap] * i, @penalties[:trailing_or_heading_sequence_gap] + (i.pred * @penalties[:consecutive_trailing_or_heading_sequence_gap])], target, i, sequence, 0, best_score)
      end
      best_alignment_obj = alignment_obj if alignment_obj && alignment_obj.first && alignment_obj.last < best_alignment_obj.last
    end
    best_alignment_obj && best_alignment_obj.last && best_alignment_obj.last <= @max_score ? best_alignment_obj : nil
  end

  def align(_target, _sequence)
    best_alignment_obj = raw_align(_target, _sequence)
    best_alignment_obj.first ? alignment_obj_to_alignment(_target, _sequence, best_alignment_obj) : nil
  end

  def alignment_obj_to_alignment(_target, _sequence, _alignment_obj)
    return nil if _alignment_obj.nil? || _alignment_obj.first.nil?
    target = _target.chars.dup
    sequence = _sequence.chars.dup
    r = _alignment_obj.first.to_s.upcase.chars.inject(['', '', '']) do |r, c|
      if c == SYMBOLS[:target_gap]
        r[0] << SYMBOLS[:gap]
        r[1] << SYMBOLS[:mismatch]
        r[2] << sequence.shift.to_s
      elsif c == SYMBOLS[:sequence_gap]
        r[0] << target.shift.to_s
        r[1] << SYMBOLS[:mismatch]
        r[2] << SYMBOLS[:gap]
      else
        r[0] << target.shift.to_s
        r[1] << c
        r[2] << sequence.shift.to_s
      end
      r
    end
    until target.empty?
      r[0] << target.shift.to_s
      r[1] << SYMBOLS[:mismatch]
      r[2] << SYMBOLS[:gap]
    end
    Alignment.new(r[0], r[1], r[2], _alignment_obj.last)
  end

  def is_gap?(s)
    return !s.nil? && (s == SYMBOLS[:trailing_or_heading_target_gap] || s == SYMBOLS[:target_gap] || s == SYMBOLS[:sequence_gap] || s == SYMBOLS[:gap])
  end

  def _align(local_alignment_obj, _target_bytes, i, _sequence_bytes, j, current_max_score)
    if local_alignment_obj.last > current_max_score
      return local_alignment_obj
    elsif j < _sequence_bytes.size
      if i >= _target_bytes.size
        return _align([local_alignment_obj.first + SYMBOLS[:target_gap], local_alignment_obj.last + @penalties[is_gap?(local_alignment_obj.first[-1]) ? :consecutive_gap : :gap]], _target_bytes, i, _sequence_bytes, j + 1, current_max_score)
      else
        if MATCHES[_sequence_bytes[j]].include?(_target_bytes[i])
          return _align([local_alignment_obj.first + SYMBOLS[:match], local_alignment_obj.last + @penalties[:match]], _target_bytes, i + 1, _sequence_bytes, j + 1, current_max_score)
        else
          # mismatch case
          mismatch = _align([local_alignment_obj.first + SYMBOLS[:mismatch], local_alignment_obj.last + @penalties[:mismatch]], _target_bytes, i + 1, _sequence_bytes, j + 1, current_max_score)

          # sequence gap
          sequence_gap = _align([local_alignment_obj.first + SYMBOLS[:sequence_gap], local_alignment_obj.last + @penalties[is_gap?(local_alignment_obj.first[-1]) ? :consecutive_gap : :gap]], _target_bytes, i + 1, _sequence_bytes, j, current_max_score)

          # target gap
          target_gap = _align([local_alignment_obj.first + SYMBOLS[:target_gap], local_alignment_obj.last + @penalties[is_gap?(local_alignment_obj.first[-1]) ? :consecutive_gap : :gap]], _target_bytes, i, _sequence_bytes, j + 1, current_max_score)

          if mismatch.last <= sequence_gap.last
            if mismatch.last <= target_gap.last
              return mismatch
            else
              return target_gap
            end
          else
            if target_gap.last <= sequence_gap.last
              return target_gap
            else
              return sequence_gap
            end
          end
        end
      end
    end
    local_alignment_obj
  end
end

class PairedEndsAlignmentAlgorithm < AlignmentAlgorithm
  SYMBOLS = {repeat_or_flank: '>'}.update(AlignmentAlgorithm::SYMBOLS).freeze

  PENALTIES = {repeat_or_flank: 0.0}.update(AlignmentAlgorithm::PENALTIES).freeze

  def initialize(max_score = 1.0, penalties = nil)
    super(max_score, penalties)
  end

  def align(_target, _paired_ends)
    # forward alignment
    i = 0
    forward_alignment_obj = [nil, @max_score + 0.1]
    while(i < (_target.size - (_paired_ends.forward.size + _paired_ends.reverse.size)))
      forward_alignment_obj_tmp = _align((i > 0 ? [SYMBOLS[:sequence_gap] * i, @penalties[:trailing_or_heading_sequence_gap] + (i.pred * @penalties[:consecutive_trailing_or_heading_sequence_gap])] : ['', 0.0]), _target.bytes, i,
      _paired_ends.forward.bytes, 0, @max_score + 0.1)
      forward_alignment_obj = forward_alignment_obj_tmp if forward_alignment_obj_tmp && forward_alignment_obj_tmp.first && forward_alignment_obj_tmp.last < forward_alignment_obj.last
      i += 1
    end

    if forward_alignment_obj.first && forward_alignment_obj.last <= @max_score
      # reverse alignment (reverse both, target and sequence)
      i = 0
      reverse_alignment_obj = [nil, @max_score + 0.1]
      while(i < _target.size - forward_alignment_obj.first.strip.size)
        reverse_alignment_obj_tmp = _align((i > 0 ? [SYMBOLS[:sequence_gap] * i, @penalties[:trailing_or_heading_sequence_gap] + (i.pred * @penalties[:consecutive_trailing_or_heading_sequence_gap])] : ['', 0.0]), _target.reverse.bytes, i,
        _paired_ends.reverse.reverse.bytes, 0, @max_score + 0.1)
        reverse_alignment_obj = reverse_alignment_obj_tmp if reverse_alignment_obj_tmp && reverse_alignment_obj_tmp.first && reverse_alignment_obj_tmp.last < reverse_alignment_obj.last
        i += 1
      end

      if forward_alignment_obj.first and reverse_alignment_obj.first
        forward_alignment_string = forward_alignment_obj.first.strip
        reverse_alignment_string = reverse_alignment_obj.first.strip.reverse
        alignment_string = forward_alignment_string
        alignment_string << (SYMBOLS[:repeat_or_flank] * (_target.size - (forward_alignment_string.size + reverse_alignment_string.size) + forward_alignment_string.count(SYMBOLS[:target_gap]) + reverse_alignment_string.count(SYMBOLS[:target_gap])))
        alignment_string << reverse_alignment_string

        return alignment_obj_to_alignment(_target, _paired_ends, [alignment_string, forward_alignment_obj.last + reverse_alignment_obj.last]) if reverse_alignment_obj.last <= @max_score
      end
    end
    nil
  end

  def alignment_obj_to_alignment(_target, _paired_ends, _alignment_obj)
    return nil if _alignment_obj.nil? || _alignment_obj.first.nil?
    target = _target.chars.dup
    sequence = _paired_ends.forward.chars.dup
    reverse = _paired_ends.reverse.chars.dup
    r = _alignment_obj.first.to_s.upcase.chars.inject(['', '', '']) do |r, c|
      if c == SYMBOLS[:target_gap]
        r[0] << SYMBOLS[:gap]
        r[1] << SYMBOLS[:mismatch]
        r[2] << sequence.shift.to_s
      elsif c == SYMBOLS[:sequence_gap]
        r[0] << target.shift.to_s
        r[1] << SYMBOLS[:mismatch]
        r[2] << SYMBOLS[:gap]
      elsif c == SYMBOLS[:repeat_or_flank]
        sequence = reverse
        r[0] << target.shift.to_s
        r[1] << SYMBOLS[:repeat_or_flank]
        r[2] << SYMBOLS[:repeat_or_flank]
      else
        r[0] << target.shift.to_s
        r[1] << c
        r[2] << sequence.shift.to_s
      end
      r
    end
    until target.empty?
      r[0] << target.shift.to_s
      r[1] << SYMBOLS[:mismatch]
      r[2] << SYMBOLS[:gap]
    end
    PairedEndsAlignment.new(r[0], r[1], r[2], _alignment_obj.last)
  end
end
