class FASTQ
  Q_VALUES = '!"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~'
  MIN_Q_VALUE = '!'.ord
  R1_R2_SPLITER = '#####'

  attr_reader :metas, :sequences, :qualities, :trim_status

  def initialize(*fastqs)
    @metas = []
    @sequences = []
    @qualities = []
    @trim_status = []
    if !fastqs.nil? && !fastqs.empty? && fastqs.any?
      @metas = fastqs.first.metas.clone
      @sequences = fastqs.first.sequences.clone
      @qualities = fastqs.first.qualities.clone
      @trim_status = fastqs.first.trim_status.clone

      if fastqs.size > 1
        fastqs[1..-1].each do |fastq|
          @metas.size.times do |i|
            @metas[i] << R1_R2_SPLITER << fastq.metas[i]
            @sequences[i] << R1_R2_SPLITER << fastq.sequences[i]
            @qualities[i] << R1_R2_SPLITER << fastq.qualities[i]
          end
        end
      end
    end
  end

  def add(meta, sequence, quality, trim_status = -1)
    @metas << meta
    @sequences << sequence
    @qualities << quality
    @trim_status << trim_status
  end

  def add_all(_metas, _sequences, _qualities, _trim_status)
    @metas = _metas.clone
    @sequences = _sequences.clone
    @qualities = _qualities.clone
    @trim_status = _trim_status.clone
  end

  def delete_if(&block)
    @metas.size.times do |i|
      if yield(@metas[i], @sequences[i], @qualities[i], @trim_status[i])
        @metas[i] = nil
        @sequences[i] = nil
        @qualities[i] = nil
        @trim_status[i] = nil
      end
    end
    @metas.compact!
    @sequences.compact!
    @qualities.compact!
    @trim_status.compact!
    true
  end

  def each(&block)
    @metas.size.times do |i|
      yield(@metas[i], @sequences[i], @qualities[i], @trim_status[i])
    end
  end

  def each_slice(l)
    r = []
    i = 0
    while i < @metas.size
      f = FASTQ.new
      f.add_all(@metas[i .. (i + l)], @sequences[i .. (i + l)], @qualities[i .. (i + l)], @trim_status[i .. (i + l)])
      r << f
      i = i + l
    end
    if block_given?
      r.each do |fastq|
        yield(fastq)
      end
    end
    r
  end

  def length
    @metas.length
  end

  def apply_q(q, replace_base = 'N')
    min_q = q.to_s[0].ord
    @qualities.each_with_index do |q, i|
      q.bytes.each_with_index do |qb, j|
        @sequences[i][j] = replace_base if qb < min_q
      end
    end
  end

  def trim_by_sequence_or_regexp!(sequence_or_regexp)
    @metas.size.times.inject(0) do |trimmed_count, i|
      if pos = @sequences[i].index(sequence_or_regexp)
        @sequences[i].slice!(pos .. -1)
        @qualities[i].slice!(pos .. -1)
        @trim_status[i] = pos
        trimmed_count + 1
      else
        trimmed_count
      end
    end
  end

  def self.percentify_q(q)
    '%.1f'%((1.0 - 10.0**((-(q.ord - 33)).to_f/10.0)) * 100.0)
  end
end
