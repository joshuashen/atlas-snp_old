#!/usr/bin/env ruby

homeDir=Dir.pwd
#$stderr.puts homeDir

#homeDir=File.dirname(File.expand_path($1))
$stderr.puts homeDir
outputDir=homeDir+'/output/'
concensusDir=homeDir+'/consensus/'
inputDir=homeDir+'/input/'
parDir=homeDir+'/par/'
refDir=homeDir+'/references/'
#refDir=homeDir+'/references/'
file_par=parDir+'refhomo.par'
fo_par=File.new(file_par,'r')
#read parameters
inputFile=''
concensusFile=''
refFile=''

fo_par.each do |line|
  cols=line.split(/:\s+/)
  
  parName=cols[0]
  p_content=cols[1]
  
  c=p_content.split(/\s+/)
  
  par_content=c[0]
  if parName=="input"
    inputFile=par_content
  end
  
  if parName=="consensus"
    concensusFile=par_content
    #$stderr.puts "$$$"
  end
  
  if parName=="reference"
    refFile=par_content
  end
end


file_output=outputDir+inputFile+'.homo.ref'
fo_output=File.new(file_output,'w')

file_con=concensusDir+concensusFile
#$stderr.print file_con
#$stderr.puts file_con
#a=`ls -l #{file_con}`
#$stderr.puts a

fo_con=File.new(file_con,'r')

#c=File.new("/data/pmrs/zwan/atlas-snp/cross_match_result/homoref-test/consensus/11_13_2007.fa.xm.ref.Cov.target_base_cov",'r')

file_snp_list=inputDir+inputFile
fo_snp_list=File.new(file_snp_list,'r')

file_ref=refDir+refFile
#file_ref='/data/pmrs/zwan/atlas-snp/human.build36.fa'
fo_ref=File.new(file_ref,'r')
ref=''
$seq={}
fo_ref.each do |line|
  if line=~ /^>(\S+)/
    ref=$1
    $seq[ref]=''
  else
    $seq[ref]<<line.chomp
  end
end
$stderr.puts "The array of sequence has been created!"

$hash_snp={}
$hash_snp=Hash.new {|h,k| h[k]=Hash.new}

#=begin
fo_snp_list.each do |line|
  cols=line.split(/\s+/)
  
  chr=cols[0]
  pos=cols[1].to_i
  refbase=cols[2]
  snpbase=cols[6]
  cov=cols[5].to_i
  refreads=cols[5].to_i-cols[10].to_i
  snpreads=cols[9].to_i
  
  $hash_snp[chr][pos]={}
  $hash_snp[chr][pos][:refbase]=refbase
  $hash_snp[chr][pos][:snpbase]=snpbase
  $hash_snp[chr][pos][:cov]=cov
  $hash_snp[chr][pos][:refreads]=refreads
  $hash_snp[chr][pos][:snpreads]=snpreads
end

$stderr.puts "The hash for snps has been created!"
fo_snp_list.close


######initial cross match##############3
$hash_cs={}
$hash_cs=Hash.new {|h,k| h[k]=Hash.new}

#######################Create hash for Crossmatch###############################3
fo_con.each do |line|
  cols=line.split(/\s+/)
  
  no=cols[2].to_i
  chr=cols[3]
  pos=cols[4].to_i
  
  $hash_cs[chr][pos]=no
end
$stderr.puts "Hash for crossmatch result has been created!"

fo_con.close

$report={}
$report=Hash.new {|h,k| h[k]=Hash.new}

$hash_cs.keys.sort.each do |ref|
  $hash_cs[ref].keys.sort.each do |pos|
    cover=$hash_cs[ref][pos].to_i
    next if $hash_snp[ref][pos]!=nil or cover<5
    
    if $report[ref][pos]==nil
      $report[ref][pos]={}
    end
    refbase=$seq[ref][pos-1,1].upcase
    gt="#{refbase}#{refbase}"
    $report[ref][pos][:gt]=gt
    $report[ref][pos][:cov]=cover
  end
end
fo_output.print "Chr\tPos\tGenotype\tCoverage\n"
$report.keys.sort.each do |ref|
  $report[ref].keys.sort.each do |pos|
    fo_output.puts "#{ref}\t#{pos}\t#{$report[ref][pos][:gt]}\t#{$report[ref][pos][:cov]}\n"
  end
end

fo_output.close
$stderr.puts "The result has been printed successfully!"
