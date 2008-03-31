## synonymous, non-synonymous, etc

# use BioRuby, could be replaced by simple Hash with a codon table.

# todo: output gene sequences with mutation in it -- for calculation of KaKs. 
 # question: how to deal with heterozygous SNPs?

#  input : reference sequence, exon coordinates, and SNP list

require 'bio'
require 'getoptlong'

$seq = {} # the DNA sequence of chromosomes
$genes = {}  # genes,  sorted by end

# codon table
$codon = {"GCT" => 'A', "GCC" => 'A', 'GCA' => 'A', "GCG" => 'A', 
  "CGT" => 'R', "CGC" => "R", "CGA" => 'R', 'CGG' => 'R', 'AGA' => 'R', 'AGG' => 'R', 
  "AAT" => 'N', "AAC" => 'N', 
  "GAT" => 'D', "GAC" => 'D', 
  "TGT" => 'C', "TGC" => 'C', 
  "CAA" => 'Q', "CAG" => 'Q', 
  "GAA" => 'E', "GAG" => 'E', 
  "GGT" => 'G', "GGC" => 'G', "GGA" => 'G', "GGG" => 'G', 
  "CAT" => 'H', "CAC" => 'H',
  "ATT" => 'I', "ATC" => 'I', "ATA" => 'I', 
  "ATG" => 'M', 
  "TTA" => 'L', "TTG" => 'L', "CTT" => 'L', "CTC" => 'L', "CTA" => 'L', "CTG" => 'L', 
  "AAA" => 'K', "AAG" => 'K', 
  "TTT" => 'F', "TTC" => 'F', 
  "CCT" => 'P', "CCC" => 'P', "CCA" => 'P', "CCG" => 'P', 
  "TCT" => 'S', "TCC" => 'S', "TCA" => 'S', "TCG" => 'S', "AGT" => 'S', "AGC" => 'S', 
  "ACT" => 'T', "ACC" => 'T', "ACA" => 'T', "ACG" => 'T', 
  "TGG" => 'W', 
  "TAT" => 'Y', "TAC" => 'Y', 
  "GTT" => 'V', "GTC" => 'V', "GTA" => 'V', "GTG" => 'V', 
  "TAG" => '*', "TGA" => '*', "TAA" => '*' }


def translate(codon)
  x = Bio::Sequence::NA.new(codon)
  a = x.translate # or a = x.translate.codes
  return a
end

def extractseq(chr)
  # extract coding DNA sequence for each gene on chr
  return if chr== '' or !$seq.key?(chr)
  
  $genes[chr].each do |gene|
    coding = ''
    gene.orfs.each do |orf|
      coding << $seq[chr][(orf[0]-1)..(orf[1]-1)]
    end
    if gene.strand == '+'
      gene.coding = coding
      gene.mutant = String.new(gene.coding)
    else
      gene.coding = coding.reverse.tr('ATGCatgc', 'TACGtacg') # reverse complementary
      gene.mutant = String.new(gene.coding)
    end
  end
  $seq[chr] = ''
  return 1
end


class Gene
  attr_accessor :chr, :orfs, :name, :nick, :coding, :strand, :exonStart, :exonEnd, :mutant, :cdsStart, :cdsEnd

  def initialize(chr,exonStart, exonEnd, orfs,strand,name,nick)
    @chr = chr
    @exonStart = exonStart
    @exonEnd = exonEnd
    @orfs = orfs.sort {|a,b| a[0] <=> b[0]}
    @strand = strand
    @name = name
    @nick = nick
    
    

    if @orfs.size > 0
      @cdsStart = @orfs[0][0] ##
      @cdsEnd = @orfs[-1][1]
    else
      @cdsStart = @exonStart
      @cdsEnd = @cdsStart
    end
  end
  
  def mutate(pos, base)
    status = "utr"
    if pos < @cdsStart or pos > @cdsEnd  or @cdsStart == @cdsEnd
      return status
    else
      shift = 0 
## !! need to consider - strand
      if @strand == '+'
        @orfs.each do |orf|
          if pos < orf[0] # intron
#          return "intron"
            status = "intron"
            return status
          elsif pos <= orf[1] # in this exon
            frame = (shift + (pos - orf[0]) ) % 3
          
            s = pos - frame - orf[0] + shift
            e  = s + 2
            aapos = s / 3 + 1
            oricodon = @coding[s..e]
            $stderr.puts "#{@name} #{@strand}  #{@coding.length}   #{s}  #{e}  #{oricodon}  #{base} #{frame} #{pos} #{shift}"
            if frame == 0
              mutcodon = base + @coding[s+1..e]
              @mutant[s,1] = base
            elsif frame == 1  # the second base
              mutcodon = oricodon[0,1] + base + oricodon[2,1]
              @mutant[s+1,1] = base
            elsif frame == 2
              mutcodon = oricodon[0..1] + base
              @mutant[s+2,1] = base
            end
            oriaa = translate(oricodon)
            mutaa = translate(mutcodon)
          
            if oriaa == mutaa # syn
              status = "coding-synonymous"
            else
              status = "coding-nonsynonymous(#{aapos}):#{oricodon}->#{mutcodon}:#{oriaa}->#{mutaa}"
            end
            return status
          else  # onto next exon
            shift += orf[1] - orf[0] + 1
          end
        end
      else # negative strand
        base = base.tr('ATGCatgc', 'TACGtacg') ## ! important
        @orfs.reverse.each do |orf|
          if pos > orf[1] # intron
            status = "intron"
            return status
          elsif pos >= orf[0] # in this exon
            frame = ((orf[1] - pos) + shift) % 3
            
            s = orf[1] - pos - frame + shift
            e = s + 2
            aapos = s / 3 + 1
            oricodon = @coding[s..e]
            $stderr.puts "#{@name} #{@strand}  #{@coding.length}   #{s}  #{e}  #{oricodon}  #{base} #{frame} #{pos} #{shift} #{orf.join(" ")}"
            if frame == 0
              mutcodon = base + @coding[s+1..e]
              @mutant[s,1] = base
              
            elsif frame == 1  # the second base 
              mutcodon = oricodon[0,1] + base + oricodon[2,1]
              @mutant[s+1,1] = base
              
            elsif frame == 2
              mutcodon = oricodon[0..1] + base
              @mutant[s+2,1] = base
              
            end
            oriaa = translate(oricodon)
            mutaa = translate(mutcodon)

            if oriaa == mutaa # syn                                              
              status = "coding-synonymous"
            else
              status = "coding-nonsynonymous(#{aapos}):#{oricodon}->#{mutcodon}:#{oriaa}->#{mutaa}"
            end
            return status
          else  # onto next exon
            shift += orf[1] - orf[0] + 1
          end
          
          
        end
      end
    end

    return status 
  end
end

opts = GetoptLong.new(
    ["--snplist", "-s", GetoptLong::REQUIRED_ARGUMENT],
    ["--reference", "-r", GetoptLong::REQUIRED_ARGUMENT],
    ["--genes", "-g", GetoptLong::REQUIRED_ARGUMENT],
    ["--output", "-o", GetoptLong::OPTIONAL_ARGUMENT],
    ["--help", "-h", GetoptLong::NO_ARGUMENT]
)

optHash = {}
opts.each do |opt, arg|
  optHash[opt] = arg
end

if optHash.key?("--help") 
  $stderr.puts "Usage: ruby __.rb -s snp_list -r reference.fa -g genes_exon_coordinates [ -o <output_prefix> ]"
  exit
end

if optHash.key?("--output")
  $prefix = optHash["--output"]
else
  $prefix = optHash["--snplist"]
end


File.new(optHash["--genes"], "r").each do |line|
  # UCSC knowngene.txt format
  ## note: sequence start at 0
  cols = line.split(/\s+/)
  name, ref, strand, txStart, txEnd, cdsStart, cdsEnd, numExons, exonStarts, exonEnds, misc = cols[0],cols[1],cols[2],cols[3].to_i+1,cols[4].to_i,cols[5].to_i+1, cols[6].to_i, cols[7].to_i, cols[8].split(','),cols[9].split(','),cols[10]
  
  orfs = [] # open reading frame
#  $stderr.puts exonStarts.class
  exonLeft = exonStarts[0].to_i + 1
  exonRight = exonEnds[-1].to_i 

  if cdsStart == cdsEnd ## no orf
    1
  else
    0.upto((exonStarts.size - 1 )) do |i|
      if cdsStart > exonEnds[i].to_i   or  cdsEnd < exonStarts[i].to_i + 1  # the whole exon is UTR
        1
      else
        s = [cdsStart,(exonStarts[i].to_i + 1) ].max
        e = [cdsEnd, (exonEnds[i].to_i ) ].min
        orfs << [s,e]
        if s > e 
          $stderr.puts "#{s} #{e} #{name}"
        end
      end
    end

  end
  #   def initialize(chr,txStart, txEnd, exons,strand,name,nick)

  gene = Gene.new(ref,exonLeft,exonRight,orfs,strand,name,misc)
  $genes[ref] = [] unless $genes.key?(ref)
  $genes[ref] << gene 
end

# sort
$bookmark = {}
$genes.each_key do |ref|
  $genes[ref].sort! {|a,b| a.exonEnd <=> b.exonEnd }
  $bookmark[ref] = 0
end

ref = ''
File.new(optHash["--reference"], "r").each do |line|
  if line=~ /^>(\S+)/
    extractseq(ref)
    ref = $1
    $seq[ref] = ''
  else
    $seq[ref] << line.chomp
  end
end
extractseq(ref)
    
ao = File.new($prefix + ".Annotation", "w")
so = File.new($prefix + ".gene_CDS.fa", "w")
mo = File.new($prefix + ".gene_CDS_mutant.fa", "w")

# assuming SNPs are sorted by base positions

$header = {'refName' => 0, 'coordinate' => 1, 'refBase' =>2, 'homopolymer' => 3, 'refEnv' => 4, 'coverage' => 5, 'SNPBase' =>6, 'numSNPReads' => 7}

# $queue = {}
File.new(optHash["--snplist"], "r").each do |line| 
  line.chomp!
  ao.print line + "\t"
  cols = line.split(/\s+/)
  if cols[0] == 'refName'  # header line
    0.upto(cols.size - 1) do |i|
      $header[cols[i]] = i
    end
  else 
    ref, pos, snpbase = cols[$header['refName']], cols[$header['coordinate']].to_i, cols[$header['SNPBase']]
    next unless $genes.key?(ref)
    
    status = 'inter-genic'

    $bookmark[ref].upto(($genes[ref].size - 1)) do |i|
      gene = $genes[ref][i]
      if gene.exonStart > pos  #  need optimization for speed
        next
      elsif gene.exonEnd < pos
        $bookmark[ref] = i
      else #  in the gene
        status = gene.mutate(pos,snpbase)
        ao.print "#{gene.name}:#{status};"
#        status = 'genic'
      end
    end

    if status == 'inter-genic'
      ao.print status
    end
  end
  ao.print "\n"
end

ao.close

$genes.each_key do |ref|
  $genes[ref].sort {|a,b| a.exonStart <=> b.exonStart}.each do |gene|
    next if gene.orfs.size < 1
    $stderr.puts "\n#{gene.name}"

    so.puts ">#{gene.name}\tori\t#{ref}\t#{gene.cdsStart}\t#{gene.cdsEnd}\t#{gene.strand}\t#{gene.nick}"
    so.puts "#{gene.coding}"
    
    mo.puts ">#{gene.name}\tmut\t#{ref}\t#{gene.cdsStart}\t#{gene.cdsEnd}\t#{gene.strand}\t#{gene.nick}"
    mo.puts "#{gene.mutant}"
  end
end

mo.close
so.close
