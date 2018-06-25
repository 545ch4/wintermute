require File.expand_path('sequence_helper.rb', File.dirname(__FILE__))
module AlignmentHelper
  include SequenceHelper

  def index_of_regular_expression(a, b, include_pattern = false)
     regexp = a.kind_of?(Regexp) ? a : b
     str = a.kind_of?(Regexp) ? b : a
     str.match(regexp) do |m|
        return m.pre_match.size + (include_pattern ? 0 : (m[0].size))
     end
     return nil
  end

  def rindex_of_regular_expression(a, b, include_pattern = false)
     regexp = Regexp.new("(#{(a.kind_of?(Regexp) ? a : b).inspect[1..-2]}).*$")
     str = a.kind_of?(Regexp) ? b : a
     str.match(regexp) do |m|
        return m.pre_match.size + (include_pattern ? (m[1].size) : 0)
     end
     return nil
  end
end
