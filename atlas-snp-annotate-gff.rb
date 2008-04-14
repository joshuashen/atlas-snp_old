## parse gff or gtf file

# description of gff format: http://www.sanger.ac.uk/Software/formats/GFF/GFF_Spec.shtml

# gtf:  http://mblab.wustl.edu/GTF2.html

## <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]

# DDB0232428      .       gene    238722  242663  .       +       .       ID=msh6;Name=msh6
# DDB0232428      Sequencing Center       mRNA    238722  245068  .       +       .       ID=DDB0202069;Parent=msh6;Name=DDB0202069;description=JC1V2_0_00120: Obtained from the Dictyostelium Genome Consortium at The Wellcome Trust Sanger Institute;translation_start=1;Dbxref=Genome V. 2.0 ID:JC1V2_0_00120,Genome V. 1.0 ID:JC1V1_0C0002_03547
# DDB0232428      Sequencing Center       CDS     238722  238799  .       +       .       Parent=DDB0202069
# DDB0232428      Sequencing Center       CDS     238890  241895  .       +       .       Parent=DDB0202069
# DDB0232428      Sequencing Center       CDS     241965  242559  .       +       .       Parent=DDB0202069
# DDB0232428      Sequencing Center       CDS     242975  243073  .       +       .       Parent=DDB0202069
# DDB0232428      Sequencing Center       CDS     243180  243580  .       +       .       Parent=DDB0202069
# DDB0232428      Sequencing Center       CDS     243662  245068  .       +       .       Parent=DDB0202069

class Gene 
  attr_accessor :chr, :orfs, :geneName, :desc, :strand, :coding,  :s, :e, :mutant, :cdsStart, :cdsEnd, :geneID
  
  def initialize(chr, geneID,geneName, s, e, strand, desc)
    @chr, @geneID,  @geneName, @s, @e, @strand, @desc = chr, geneID,  geneName, s, e, strand, desc
    @orfs = []
    @cdsStart, @cdsEnd = e, s  # initialize
    @coding = ''
    @mutant = ''
  end
    
  def addOrf(s,e)
    @orfs << [s,e]
    @cdsStart = s if s < @cdsStart
    @cdsEnd = e if e > @cdsEnd
    @orfs.sort! {|a,b| a[0] <=> b[0]}
  end
  
  def mutate(pos, base)
    status = "utr"
    if pos < @cdsStart or pos > @cdsEnd  or @cdsStart == @cdsEnd
      return status
    else
      shift = 0
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
            #              $stderr.puts "#{@name} #{@strand}  #{@coding.length}   #{s}  #{e}  #{oricodon}  #{base} #{frame} #{pos} #{shift}"
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
          #  $stderr.puts "#{@geneName} #{@strand}  #{@coding.length}   #{s}  #{e}  #{oricodon}  #{base} #{frame} #{pos} #{shift} #{orf.join(" ")}"
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
  end # end of mutate()
end  # end of class Gene

class Gff
  attr_accessor :genes
  def initialize(file)
    @file = file

    @genes = {}
  end
  
  def process()

    gene = nil
    desc = ''
    flag = 0
    File.new(@file,'r').each do |line|
      next if line=~ /^#/
        cols=line.split(/\t/)  # fields are separated by tabs
      next if cols.size < 9 # at least 9 columns
      
      ref, source, feature, s, e, strand, frame, info = cols[0], cols[1], cols[2], cols[3].to_i, cols[4].to_i, cols[6], cols[7], cols[8]
      
      attr = {}
      info.split(';').each do |tt|

        x =tt.split('=')
#        $stderr.puts x[0]
        attr[x[0]] = x[1]
#        $stderr.puts attr['ID']
      end
      
      if feature == "gene"  # a new gene
        desc = attr["description"]
      elsif feature == "mRNA" # a new gene product
        geneID= attr["ID"]
        parent = attr["Parent"]
        gene = Gene.new(ref, geneID, parent, s, e, strand,desc)
 #       $stderr.puts "#{ref}\t#{geneID}"
        @genes[ref] = [] if !@genes.key?(ref)
        @genes[ref] << gene          
#        $stderr.puts "#{ref}\t#{geneID}"
      elsif feature == "CDS"
        gene.addOrf(s,e)
      end
    end
#    $stderr.puts @genes.size
    return 1
  end
end  # end of class Gff


require 'bio'
require 'getoptlong'
        
$seq = {} # the DNA sequence of chromosomes    
$genes = nil  # genes,  sorted by end

def translate(codon)
  x = Bio::Sequence::NA.new(codon)
  a = x.translate # or a = x.translate.codes 
  return a
end

def extractseq(chr)
  # extract coding DNA sequence for each gene on chr
  return if chr== '' or !$genes.key?(chr) or !$seq.key?(chr)
#  $stderr.puts "#{chr}"
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

gff = Gff.new(optHash["--genes"])
gff.process
$genes = gff.genes

# $stderr.puts "genes:  #{$genes.keys.join(' ')}"
$bookmark = {}
$genes.each_key do |ref|
  $genes[ref].sort! {|a,b| a.e <=> b.e }
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

$header = {'refName' => 0, 'coordinate' => 1, 'refBase' =>2, 'homopolymer' => 3, 'refEnv' => 4, 'coverage' => 5, 'SNPBase'  =>6, 'numSNPReads' => 7}

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
    
    if !$genes.key?(ref)
      status = 'inter-genic'
    else
    
      status = 'inter-genic'

      $bookmark[ref].upto(($genes[ref].size - 1)) do |i|
        gene = $genes[ref][i]
        if gene.s > pos  #  need optimization for speed                                                              
          next
        elsif gene.e < pos
          $bookmark[ref] = i
        else #  in the gene                                                                                                  
          status = gene.mutate(pos,snpbase)
          ao.print "#{gene.geneName}:#{status};"
          #        status = 'genic'                                                                                                  
        end
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
  $genes[ref].sort {|a,b| a.s <=> b.s}.each do |gene|
    next if gene.orfs.size < 1
#    $stderr.puts "\n#{gene.geneID}"

    so.puts ">#{gene.geneID}\tori\t#{ref}\t#{gene.cdsStart}\t#{gene.cdsEnd}\t#{gene.strand}\t#{gene.geneName}"
    ss = gene.coding.size
    0.upto(ss/80) do |i|
      so.print(gene.coding[i*80..(i+1)*80-1] + "\n")
    end
    
      
#    so.puts "#{gene.coding}"

    mo.puts ">#{gene.geneID}\tmut\t#{ref}\t#{gene.cdsStart}\t#{gene.cdsEnd}\t#{gene.strand}\t#{gene.geneName}"
    ss = gene.mutant.size
     0.upto(ss/80) do |i|
      mo.print(gene.mutant[i*80..(i+1)*80-1] + "\n")
    end
 
#    mo.puts "#{gene.mutant}"
  end
end

mo.close
so.close
