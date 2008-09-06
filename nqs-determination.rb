#!/usr/bin/env ruby

##Intialize directory###############
package=File.dirname(File.expand_path($0))
$stderr.puts package
homeDir=Dir.pwd

inputDir=homeDir+'/input/'
parDir=homeDir+'/par/'

outputDir=homeDir+'/output/'

##Intialize files##################
parFile=parDir+'nqs.par'
fo_par=File.new(parFile,'r')

prefixName=''
fo_par.each do |line|
  cols=line.split(/:\s+/)
  
  parName=cols[0]
  pContent=cols[1]
  
  c=pContent.split(/\s+/)
  parContent=c[0]
  
  if parName=="qual"
    prefixName=parContent
  end
end

file_dis_nqs=outputDir+prefixName+'.nqs.dis'
fo_dis_nqs=File.new(file_dis_nqs,'w')


###Start to do nqs for each batch###############################
str="ruby #{package}/nqs-pass-distribution.rb -i $f"
cmd="count=0; for f in #{inputDir}#{prefixName}_*.fa.qual; do echo $f; echo $count; count=`expr $count + 1`; bsub -J #{prefixName}_$count -o #{outputDir}#{prefixName}_$count.o -e #{outputDir}#{prefixName}_$count.e #{str}; done"
#$stderr.puts cmd

system(cmd)
temp=outputDir+'temp'
fo_temp=File.new(temp,'w')
fo_temp.print "Read_name\tNo_passed\tNo_total\tPassed_ratio\n"
fo_temp.close

##Check the job status##############################
loop {
  bjobsOutput=`bjobs`
  bj=(bjobsOutput !~ /#{prefixName}_/)
  #$stderr.puts bj
  if (bjobsOutput.empty? or bjobsOutput !~ /#{prefixName}_/ )
    
    #####################################
    contains=Dir.new("output").entries
    #$stderr.puts contains
    b=[]
    failed=false
    while contains.size>0
      a=contains.shift
      
      if a=~ /#{prefixName}_[0-9]\.o/ or a=~ /#{prefixName}_[0-9][0-9]\.o/
         b<<a 
      end
    end
    #$stderr.puts a.size
    while b.size>0
      out=b.shift
      grepResult=`grep "Exit" #{outputDir}#{out}`
      if grepResult!=""
        if out=~ /^(\S+)\.(\S+)/
          name=$1
          type=$2
          $stderr.puts name+" failed! Please check "+name+".e for error information!"
          failed=true
          exit
        end 
      end
    end
    
    ######################################
    if !failed
      `cat #{temp} #{outputDir}#{prefixName}_*.fa.qual.nqs > #{outputDir}#{prefixName}.nqs`
      $stderr.puts "successful!"
    end
    break
  else
    #$stderr.puts "wait!!!!"
  end
  sleep(60)
}
$stderr.puts "nqs result has been created!"
############################################
cols=[]
count=[]
0.upto(39) do |i|
  count<<0
end

dis_list=`ls -d #{outputDir}#{prefixName}_*.fa.qual.nqs.dis`
cols=dis_list.split(/\n/)

0.upto(cols.size-1) do |i|
  name=cols[i]
  n=0
  File.new(name,'r').each do |line|
    c=line.split(/\s+/)
    no=c[1].to_i
    count[n]+=no
    n+=1
  end
end

fo_dis_nqs.puts "PassRatio_range\tNo_of_reads"

#c=""
0.upto(9) do |i|
  0.upto(3) do |j|
  
  a1=count.shift
  j1=j+1
  fo_dis_nqs.print "#{i.to_f/10+j.to_f/40}-#{i.to_f/10+j1.to_f/40}\t#{a1}\n"
  end
end
fo_dis_nqs.close

$stderr.puts "Distribution has been created!"

##Start to delete temp files#####################
`rm #{outputDir}#{prefixName}_*.fa.qual.nqs.dis`
`rm #{outputDir}#{prefixName}_*.fa.qual.nqs`
`rm #{outputDir}#{prefixName}_*.o`
`rm #{outputDir}#{prefixName}_*.e`
`rm #{temp}`