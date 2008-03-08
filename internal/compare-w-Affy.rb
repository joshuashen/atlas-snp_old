## compare SNP calls with Affy data

class Affy
  attr_accessor :ref, :pos, :genotype, :alleles, :rs
  
  def initialize(ref,pos,rs, gtype,alleles)
    @ref, @pos, @genotype, @alleles,@rs = ref, pos, gtype, alleles,rs
  end
end


require 'getoptlong'

opts = GetoptLong.new(
    ["--affy", "-a" , GetoptLong::REQUIRED_ARGUMENT],
    ["--snp", "-s", GetoptLong::REQUIRED_ARGUMENT],
    ["--help", "-h", GetoptLong::NO_ARGUMENT]
)

optHash={}
opts.each do |opt,arg|
  optHash[opt] = arg
end

if optHash.key?("--help") or !optHash.key?("--affy") or !optHash.key?("--snp")
  $stderr.puts "Usage: ruby __.rb -s snp.list -a affy.genotype > output"
  exit
end


$genotypes =Hash.new {|h,k| h[k]={} }


File.new(optHash["--affy"], "r").each do |line|
  cols = line.split(/\s+/)
  ref,pos,rs, gtype = cols[1],cols[2].to_i, cols[3], cols[-1]
  
  if gtype =~ /^SAME/
    genotype = 'ref'
    alleles = ''
  elsif gtype =~ /^SNPHOM\:(\S)\>(\S)/
    genotype = 'hom'
    alleles = $2
  elsif gtype =~ /^SNPHET\:(\S)\>(\S)\,(\S)/
    genotype = "het"
    alleles = $2 + '_' + $3
  end

  affy = Affy.new(ref,pos,rs, genotype, alleles)
  $genotypes[ref][pos] = affy
end

File.new(optHash["--snp"], "r").each do |line|
  cols = line.split(/\s+/)
  ref, pos, snpbase, cov, num, pr  = cols[0],cols[1].to_i,cols[6], cols[5].to_i, cols[9].to_i, cols[13].to_f

  if $genotypes.key?(ref) and $genotypes[ref].key?(pos)
    puts line.chomp! + "\taffy:#{$genotypes[ref][pos].genotype}\t#{$genotypes[ref][pos].alleles}"
  else
    puts line.chomp! + "\taffy:none\tNA"
  end
end
