## convert encode info to a format that can be used for atlas-mapper-coverage.rb

# original format:
# ENm010  7       200000  27124046        27224045        BCM     Old     HOXA_cluster

# target format:
#  target_name reference_name start_on_ref end_on_ref direction(1 or -1)

while line=ARGF.gets do 
  cols = line.split(/\s+/)
  encode, chr, l, s, e, info = cols[0],cols[1],cols[2],cols[3],cols[4],cols[5..-1]
  
  puts "#{encode}\tchr#{chr}\t#{s}\t#{e}\t1\t#{info.join(',')}"
end

