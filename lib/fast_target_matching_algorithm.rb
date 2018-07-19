require File.expand_path('target.rb', File.dirname(__FILE__))
require File.expand_path('alignment.rb', File.dirname(__FILE__))
require 'levenshtein'

class FastTargetMatchingAlgorithm < AlignmentAlgorithm
  DIFF_AROUND_EXACT_MATCH = 1
  DECISIONS = {
    :unkown_reason => false,
    :match => true,
    :no_forward_partial_exact_match => false,
    :forward_partial_exact_match_not_within_range => false,
    :forward_partial_exact_match_too_many_Ns => false,
#    :r2_detected_but_no_reverse_partial_exact_match => false,
    :no_reverse_partial_exact_match => false,
    :reverse_partial_exact_match_too_many_Ns => false,
    :forward_partial_exact_match_after_reverse_partial_exact_match => false,
    :no_forward_alignment => false,
    :no_reverse_alignment => false,
    :no_forward_partial_exact_match_due_to_second_half_mismatch => false,
    :no_forward_partial_exact_match_due_to_first_half_mismatch => false
  }

  def initialize(_max_score = 1.0, _penalties = nil, _options = {})
    super(_max_score, _penalties)
    @options = {:max_n => 3, :match_forward => true, :match_reverse => true, :primer_trimming => true, :max_distance_of_forward_sequence_from_start => -1, :max_primer_half_levenshtein_distance => 1}.update(_options)
  end

  def match(target, sequence)
    (reason, forward_pos, forward_alignment, reverse_pos, reverse_alignment) = raw_match(target, sequence)
    @options[:primer_trimming] ? [reason, (DECISIONS[reason] ? sequence[(forward_pos + forward_alignment.gsub('T', '').size) ... reverse_pos] : nil)] : [reason, (DECISIONS[reason] ? sequence[forward_pos ... (reverse_pos + reverse_alignment.gsub('T', '').size)] : nil)]
  end

  def match_with_alignment(target, sequence)
    (reason, forward_pos, forward_alignment, reverse_pos, reverse_alignment) = raw_match(target, sequence)
    return [reason, nil] unless DECISIONS[reason]

    alignment_primers = ''
    forward_primer_copy = target.forward_primer.original_sequence.dup
    reverse_primer_copy = target.reverse_primer.original_sequence.dup
    alignment = ''
    alignment_sequence = ''
    sequence_copy = sequence.dup

    forward_alignment.each_char do |alignment_char|
      if alignment_char == SYMBOLS[:sequence_gap]
        alignment_primers += SYMBOLS[:gap]
        alignment_sequence += sequence_copy.slice!(0)
      elsif alignment_char == SYMBOLS[:target_gap]
        alignment_primers += forward_primer_copy.slice!(0)
        alignment_sequence += SYMBOLS[:gap]
      elsif alignment_char == SYMBOLS[:match]
        alignment_primers += forward_primer_copy.slice!(0).downcase
        alignment_sequence += sequence_copy.slice!(0).downcase
      else
        alignment_primers += forward_primer_copy.slice!(0)
        alignment_sequence += sequence_copy.slice!(0)
      end
      alignment += alignment_char
    end
    
    alignment_primers += '@'
    alignment += '@'
    alignment_sequence += '@'

    ((forward_pos + forward_alignment.gsub('T', '').length) ... reverse_pos).each do |i|
      alignment_primers += SYMBOLS[:gap]
      alignment += SYMBOLS[:gap]
      alignment_sequence += sequence_copy.slice!(0)
    end

    alignment_primers += '@'
    alignment += '@'
    alignment_sequence += '@'

    reverse_alignment.each_char do |alignment_char|
      if alignment_char == SYMBOLS[:sequence_gap]
        alignment_primers += SYMBOLS[:gap]
        alignment_sequence += sequence_copy.slice!(0)
      elsif alignment_char == SYMBOLS[:target_gap]
        alignment_primers += reverse_primer_copy.slice!(0)
        alignment_sequence += SYMBOLS[:gap]
      elsif alignment_char == SYMBOLS[:match]
        alignment_primers += reverse_primer_copy.slice!(0).downcase
        alignment_sequence += sequence_copy.slice!(0).downcase
      else
        alignment_primers += reverse_primer_copy.slice!(0)
        alignment_sequence += sequence_copy.slice!(0)
      end
      alignment += alignment_char
    end
#    puts alignment_primers
#    puts alignment
#    puts alignment_sequence

    return_repeat = sequence[(forward_pos + forward_alignment.gsub('T', '').size) ... reverse_pos]
    
    return [reason, return_repeat, alignment_primers, alignment, alignment_sequence]
  end

private

  def raw_match(target, sequence)
    sequence_bytes = sequence.bytes

    # idea: one part of primer should match perfectly
    # todo: now there are only halfs found to match perfectly => split into mor parts if neccessary e.g. long primers
    found_forward_pos = nil
    found_reverse_pos = nil
    pos = nil

    no_forward_match_reason = nil
    if target.empty_forward_primer? || !@options[:match_forward]
      found_forward_pos = 0
    elsif pos = sequence.index(target.forward_primer.sequence)
      found_forward_pos = pos
    elsif pos = sequence.index(target.forward_primer.first_half)
      if Levenshtein.distance(sequence[(pos + target.forward_primer.first_half_length)...(pos + target.forward_primer.length)], target.forward_primer.original_second_half) <= @options[:max_primer_half_levenshtein_distance]
        found_forward_pos = pos
      else
        no_forward_match_reason = [:no_forward_partial_exact_match_due_to_second_half_mismatch]
      end
    elsif pos = sequence.index(target.forward_primer.second_half)
      if Levenshtein.distance(sequence[(pos - target.forward_primer.first_half_length)...pos], target.forward_primer.original_first_half) <= @options[:max_primer_half_levenshtein_distance]
        found_forward_pos = (pos + target.forward_primer.second_half_length) - target.forward_primer.length
      else
        no_forward_match_reason = [:no_forward_partial_exact_match_due_to_first_half_mismatch]
      end
    end

    if found_forward_pos.nil?
      return no_forward_match_reason || [:no_forward_partial_exact_match]
    else
      if @options[:max_n] && sequence[found_forward_pos...(found_forward_pos + target.forward_primer.length)].count('N') > @options[:max_n]
        return [:forward_partial_exact_match_too_many_Ns]
      elsif @options[:max_distance_of_forward_sequence_from_start] >= 0 && found_forward_pos > @options[:max_distance_of_forward_sequence_from_start]
        return [:forward_partial_exact_match_not_within_range]
      else
        # check r2 match if availabe and desired
        # if @options[:scan_for_r2] && (pos = sequence.index(FASTQ::R1_R2_SPLITER))
        #   r2 = sequence[(pos + FASTQ::R1_R2_SPLITER.length) .. -1]
        #   if r2.index(target.reverse_primer.first_half).nil? && r2.index(target.reverse_primer.second_half).nil?
        #     return [:r2_detected_but_no_reverse_partial_exact_match]
        #   end
        # end

        if target.empty_reverse_primer? || !@options[:match_reverse]
          found_reverse_pos = sequence.length
        elsif pos = sequence.index(target.reverse_primer.sequence)
          found_reverse_pos = pos
        elsif pos = sequence.index(target.reverse_primer.first_half)
          found_reverse_pos = pos
        elsif pos = sequence.index(target.reverse_primer.second_half)
          found_reverse_pos = (pos - target.reverse_primer.length) + target.reverse_primer.first_half_length
        end

        # we found a potential match! check with real alignemnt
        if found_reverse_pos.nil?
          return [:no_reverse_partial_exact_match]
        elsif @options[:max_n] && sequence[found_reverse_pos...(found_reverse_pos + target.reverse_primer.length)].count('N') > @options[:max_n]
          return [:reverse_partial_exact_match_too_many_Ns]
        else
          if found_forward_pos >= found_reverse_pos
            return [:forward_partial_exact_match_after_reverse_partial_exact_match]
          else
            if target.empty_forward_primer? || !@options[:match_forward]
              best_forward_pos = found_forward_pos
            else
              best_forward_score = @max_score + 0.1
              best_forward_alignment_obj = [nil, best_forward_score]
              best_forward_pos = found_forward_pos
              ((found_forward_pos - DIFF_AROUND_EXACT_MATCH) .. (found_forward_pos + DIFF_AROUND_EXACT_MATCH)).each do |i|
                forward_alignment_obj = _align(['', 0.0], sequence_bytes, i, target.forward_primer.original_sequence_as_bytes, 0, best_forward_score)
                if forward_alignment_obj && forward_alignment_obj.first && forward_alignment_obj.last < best_forward_alignment_obj.last
                  best_forward_alignment_obj = forward_alignment_obj
                  best_forward_pos = i
                end
              end
            end

            # verified partial exact macth with real alignment match (forward) => let's check the reverse
            if @options[:match_forward] && !target.empty_forward_primer? && (best_forward_alignment_obj.nil? || best_forward_alignment_obj.last.nil? || best_forward_alignment_obj.last > @max_score)
              return [:no_forward_alignment]
            else
              if target.empty_reverse_primer? || !@options[:match_reverse]
                best_reverse_pos = found_reverse_pos
              else
                best_reverse_score = @max_score + 0.1
                best_reverse_alignment_obj = [nil, best_reverse_score]
                best_reverse_pos = found_reverse_pos
                ((found_reverse_pos - DIFF_AROUND_EXACT_MATCH) .. (found_reverse_pos + DIFF_AROUND_EXACT_MATCH)).each do |i|
                  reverse_alignment_obj = _align(['', 0.0], sequence_bytes, i, target.reverse_primer.original_sequence_as_bytes, 0, best_reverse_score)
                  if reverse_alignment_obj && reverse_alignment_obj.first && reverse_alignment_obj.last < best_reverse_alignment_obj.last
                    best_reverse_alignment_obj = reverse_alignment_obj
                    best_reverse_pos = i
                  end
                end
              end

              # found both alignemnts
              if @options[:match_reverse] && !target.empty_reverse_primer? && (best_reverse_alignment_obj.nil? || best_reverse_alignment_obj.last.nil? || best_reverse_alignment_obj.last > @max_score)
                return [:no_reverse_alignment]
              else
                return [:match, best_forward_pos, ((!@options[:match_forward] || target.empty_forward_primer?) ? '' : best_forward_alignment_obj.first), best_reverse_pos, ((!@options[:match_reverse] || target.empty_reverse_primer?) ? '' : best_reverse_alignment_obj.first)]
              end
            end
          end
        end
      end
    end
    return [:unkown_reason]
  end
end
