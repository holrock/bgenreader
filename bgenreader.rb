require 'zlib'
require 'stringio'

class BgenReader
  def self.open(filename, &block)
    BgenReader.new(filename).each &block
  end

  def initialize(filename)
    @filename = filename
  end

  def each
    File.open(@filename, "rb:ASCII-8BIT") do |f|
      buf = "\x0" * 4096
      @header = parse_header_block(f, buf)
      # skip sample data
      f.seek(4 + @header[:variant_block_offset], IO::SEEK_SET)

      header[:num_variant_blocks].times do
        vid = parse_variant_id_data(f, buf)
        data_length = read_uint(f, buf)
        uncompress_data_length = read_uint(f, buf)
        data_length -= 4 unless @header[:compressed_snp_blocks] == 0
        f.read(data_length, buf)
        data = Zlib::Inflate.inflate(buf)
        sio = StringIO.new(data)
        prob = parse_genotype_data_layout2(sio, buf)
        yield(vid, prob)
      end
    end
  end

  private

  USHORT_MAX = 2 ** 16
  UINT_MAX = 2 ** 32

  def read_uint(f, buf = nil)
    f.read(4, buf)
    val = buf.unpack1("V")
    raise "parse error: overflow uint #{val}" if val > UINT_MAX
    val
  end

  def read_ushort(f, buf)
    f.read(2, buf)
    val = buf.unpack1("v")
    raise "parse error: overflow ushort #{val}" if val > USHORT_MAX
    val
  end

  def parse_header_block(f, buf)
    variant_block_offset = read_uint(f, buf)
    header_block_length = read_uint(f, buf)
    num_variant_blocks = read_uint(f, buf)
    num_samples = read_uint(f, buf)
    magic = f.read(4)
    f.read(header_block_length - 20, buf) # free data area
    flags = read_uint(f, buf)
    compressed_snp_blocks = flags & 0b11
    layout = (flags >> 2) & 0b1111
    sample_identifiers = (flags >> 31) & 0b1

    raise "layout #{layout} not supported yet" unless layout == 2

    {
      variant_block_offset: variant_block_offset,
      header_block_length: header_block_length,
      num_variant_blocks: num_variant_blocks,
      num_samples: num_samples,
      magic: magic,
      flags: flags.to_s(2),
      compressed_snp_blocks: compressed_snp_blocks,
      layout: layout,
      sample_identifiers: sample_identifiers,
    }
  end

  def parse_variant_id_data(f, buf)
    variant_id_length = read_ushort(f, buf)
    variant_id = f.read(variant_id_length)

    rsid_length = read_ushort(f, buf)
    rsid = f.read(rsid_length)

    chr_length = read_ushort(f, buf)
    chr = f.read(chr_length)

    variant_pos = read_uint(f, buf)
    num_alleles = read_ushort(f, buf)
    alleles = num_alleles.times.map do
      len = read_uint(f, buf)
      f.read(len)
    end

    {
      variant_id_length: variant_id_length,
      rsid_length: rsid_length,
      chr_length: chr_length,
      variant_pos: variant_pos,
      num_alleles: num_alleles,
      variant_id: variant_id,
      rsid: rsid,
      chr: chr,
      alleles: alleles
    }
  end

  def parse_genotype_data_layout2(f, buf)
    num_individuals = read_uint(f, buf)
    num_alleles = read_ushort(f, buf)
    min_ploidy, max_ploidy, = f.read(2, buf).bytes
    ploidy_list = f.read(num_individuals)
    phased = f.read(1).bytes.first
    num_bits = f.read(1).bytes.first
    raise "#{num_bits}bits data not supported yet" unless num_bits == 8

    d = ((2 ** num_bits) - 1).to_f
    case phased
    when 0
      prob = ploidy_list.each_byte.map do |ploidy|
        next nil if ploidy > 63
        al = f.read(ploidy * (num_alleles - 1)).each_byte.map { |g| g / d }
        al << 1.0 - al.sum
      end
    when 1
      raise "phased data not supported yet"
    else
      raise "Unkown phase flag #{phased}"
    end

    {
      num_individuals: num_individuals,
      num_alleles: num_alleles,
      min_ploidy: min_ploidy,
      max_ploidy: max_ploidy,
      phased: phased,
      num_bits: num_bits,
      prob_size: prob.size
    }
  end
end