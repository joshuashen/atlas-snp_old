## model: 
##  the duplicated regions must have:
##   1. significantly higher coverage than average
##   2. increase of "heterozygosity", i.e., the ratio of SNP reads / coverage is lower; 
##   3. heterogenity of reads: the reads show SNPs have overall higher mismatch rate
##   4. tails

class Snp
  attr_accessor :chr, :pos, :refbase, :cov, :snpbase, :num, :pr

  def initialize(chr,pos,refbase,cov,snpbase,num,pr)
    @chr,@pos,@refbase,@cov,@snpbase, @num, @pr = chr,pos.to_i,refbase,cov.to_i,snpbase,num.to_i,pr.to_f
    
  end
end

class Gene
  attr_accessor :name, :snps, :coverage

  def initialize()
    @name = ''
    @snps = {}
    @coverage = []
  end

  def slide
# collect the snp information in sliding windows
# the window size: 2k;  overlapping 1k 
# information to collect:  SNP rate, het_SNP rate, and average coverage
    genesize = @coverage.size
   #  $stderr.puts "#{genesize}\t#{@snps.size}"
    steps = (genesize / 1000.0).floor
    str = ''
    0.upto(steps-1) do |i|
      window1 = Range.new(i*1000 + 1, [genesize, (i+1)*1000].min)
      window2 = Range.new(i*1000 + 501, [genesize, (i+1)*1000 + 500].min)
      str << "#{@name}\t" + collectInfo(window1) + "\n"
      str << "#{@name}\t" + collectInfo(window2) + "\n"
    end
    lastwindow = Range.new(steps*1000 + 1, [steps*1000+500, genesize].min)
    str << "#{@name}\t" + collectInfo(lastwindow) + "\n"
    return str
  end

  def collectInfo(range)
    s,e =range.first,range.last
   # $stderr.puts "#{s}\t#{e}"
    nonzero = @coverage[s-1..e-1].select {|i| i>0}
    t = 0
    nonzero.map {|i| t+= i}
    avgCov = t / ( e -s + 1.0)
    snps = @snps.select {|k,v| k>=s and k<=e}
    hetSNP = snps.select {|i| i[1].num < i[1].cov*0.6}.size
    return "#{s}\t#{e}\t#{avgCov}\t#{snps.size/(nonzero.size+0.001)}\t#{hetSNP/(snps.size+0.001)}"
  end
  
end

snpf = ARGV[0]
covf = ARGV[1]


$snps = Hash.new {|h,k| h[k]={}}

File.new(snpf,'r').each do |line|
  cols = line.split(/\s+/)
  if cols[0] != 'refName'
    snp = Snp.new(cols[0], cols[1], cols[2],cols[5],cols[6],cols[9],cols[-1])
    $snps[cols[0]][cols[1].to_i] = snp
  end
end

# $stderr.puts $snps.size 

gene = Gene.new()


File.new(covf, 'r').each do |line|
  cols = line.split(/\s+/)
  name,pos, cov, ref, cor = cols[0], cols[1].to_i, cols[2].to_i, cols[3], cols[4].to_i
#  $stderr.puts "#{ref}\t#{cor}" 
  if name != gene.name # meet a new gene
    puts gene.slide if gene.name != ''
    gene = Gene.new()
    gene.name = name
  end

  gene.coverage[pos-1] = cov
  if $snps.key?(ref) and $snps[ref].key?(cor) # a snp position

    gene.snps[pos] = $snps[ref][cor]
  end
end    
puts gene.slide

exit
