refList = ARGV[0]

File.new(refList, 'r').each do |line|
  if line=~ /^(\S+)/
    ref =  $1
    cmd = "time ruby ~/work/ruby/atlas-snp/atlas-mapper.rb -r #{ref} -q div_9.fa -z"
    system(cmd)
  end
end


