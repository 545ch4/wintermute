class FASTQ
  Q_VALUES = '!"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~'
  MIN_Q = '!'
  MIN_Q_VALUE = MIN_Q.ord
  MINIMAL_MERGE_OVERLAP = 30
  MISSING_OVERLAP_SEPERATOR_SEQUENCE = 'N'*10
  MISSING_OVERLAP_SEPERATOR_QUALITY = '~'*10

  attr_reader :metas, :sequences, :qualities

  def initialize(*fastqs)
    @metas = []
    @sequences = []
    @qualities = []
    if !fastqs.nil? && !fastqs.empty? && fastqs.any?
      @metas = fastqs.first.metas.clone
      @sequences = fastqs.first.sequences.clone
      @qualities = fastqs.first.qualities.clone

      if fastqs.size > 1
        fastqs[1..-1].each do |fastq|
          @metas.size.times do |i|
            @metas[i] << fastq.metas[i]
            @sequences[i] << fastq.sequences[i]
            @qualities[i] << fastq.qualities[i]
          end
        end
      end
    end
  end

  def add(meta, sequence, quality)
    @metas << meta
    @sequences << sequence
    @qualities << quality
  end

  def add_all(_metas, _sequences, _qualities)
    @metas = _metas.clone
    @sequences = _sequences.clone
    @qualities = _qualities.clone
  end

  def delete_if(&block)
    @metas.size.times do |i|
      if yield(@metas[i], @sequences[i], @qualities[i])
        @metas[i] = nil
        @sequences[i] = nil
        @qualities[i] = nil
      end
    end
    @metas.compact!
    @sequences.compact!
    @qualities.compact!
    true
  end

  def each(&block)
    @metas.size.times do |i|
      yield(@metas[i], @sequences[i], @qualities[i])
    end
  end

  def each_slice(l)
    r = []
    i = 0
    while i < @metas.size
      f = FASTQ.new
      f.add_all(@metas[i .. (i + l)], @sequences[i .. (i + l)], @qualities[i .. (i + l)])
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
        trimmed_count + 1
      else
        trimmed_count
      end
    end
  end

  def self.percentify_q(q)
    '%.1f'%((1.0 - 10.0**((-(q.ord - 33)).to_f/10.0)) * 100.0)
  end

  def merge!(fastq)
    return false if @metas.size != fastq.metas.size
    @metas.size.times do |i|
      merged = false
      if @sequences[i].length >= MINIMAL_MERGE_OVERLAP || fastq.sequences[i].length >= MINIMAL_MERGE_OVERLAP
        fastq_short = @sequences[i].length > fastq.sequences[i].length ? fastq : self
        fastq_long = @sequences[i].length > fastq.sequences[i].length ? self : fastq
        (MINIMAL_MERGE_OVERLAP .. fastq_short.sequences[i].length).to_a.reverse.each do |overlapp_length|
          if pos = fastq_long.sequences[i].index(fastq_short.sequences[i][0 ... overlapp_length])
            @sequences[i] = (fastq_long.sequences[i][0 ... pos] << fastq_short.sequences[i])
            new_quality = (fastq_long.qualities[i][0 ... pos] << fastq_short.qualities[i])
            (0 ... overlapp_length).each do |j|
#              puts "new_quality.length => #{new_quality.length}, fastq_short.qualities[i].length => #{fastq_short.qualities[i].length}, overlapp_length => #{overlapp_length}, j => #{j}, pos => #{pos}"
              new_quality[pos + j] = fastq_short.qualities[i][j] if new_quality[pos + j].ord < fastq_short.qualities[i][j].ord
            end
            @qualities[i] = new_quality
            @metas[i] << " merged successfully with " << fastq.metas[i]
            merged = true
            break
          end
        end
      end

      unless merged
        if @sequences[i].length >= fastq.sequences[i].length
          @metas[i] << ' untouced due to size mismatch'
        else
          @metas[i] << ' replaced due to size mismatch with ' << fastq.metas[i]
          @sequences[i] = fastq.sequences[i]
          @qualities[i] = fastq.qualities[i]
        end
      end
    end
  end
end
