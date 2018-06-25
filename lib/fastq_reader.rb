require 'stringio'
require 'zlib'
require File.expand_path('fastq.rb', File.dirname(__FILE__))
require File.expand_path('sequence_helper.rb', File.dirname(__FILE__))

class FASTQReader
  def self.read(filename, options = {})
    fastq = FASTQ.new
    File.open(filename) do |raw_io|
      loop do
	if File.extname(filename) != '.gz'
	    io = raw_io
	else
    	    io = Zlib::GzipReader.new raw_io
	end
        while !io.eof? do
          begin
            meta = io.readline.strip
            sequence = io.readline.strip
            plus = io.readline.strip[0]
            quality = io.readline.strip

            if options[:reverse_complement]
              sequence = SequenceHelper::reverse_complement(sequence)
            end

            if plus == '+'
              fastq.add(meta, sequence, quality)
            else
              puts "meta = #{meta.inspect}"
              puts "sequence = #{sequence.inspect}"
              puts "plus = #{plus.inspect}"
              puts "quality = #{quality.inspect}"
              exit 1
            end
          rescue EOFError
            print " unexpected EOF in #{filename} "
          end
        end
	if io.kind_of?(Zlib::GzipReader)
    	    unused = io.unused
    	    io.finish
    	    break if unused.nil?
            raw_io.pos -= unused.length
	else
	    break
	end
      end
    end
    fastq
  end
end
