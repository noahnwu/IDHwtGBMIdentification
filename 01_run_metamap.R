################### CHECK THAT THE DIRECTORIES NEEDED FOR ANALYSIS EXIST ################### 
#location of metamap binary bin
ini = readRDS("ini.RDS")
metamap.loc = ini$metamap.loc
#location of input notes (in the directory already)
input.dir = ini$input.dir
#location for metamap to leave output
output.dir = ini$metmap.out


################### GET VECTOR OF INPUT FILES FOR METAMAP ################### 
print("Identifying Input Files...")
input.files = dir(input.dir, full.names = T)
output.files = paste0(output.dir, dir(input.dir, full.names = F ), ".out")

################### SPECIFY THE BASE COMMAND TO BE RUN ON EACH FILE ################### 
metamap.base = paste0(metamap.loc, "metamap --silent -y --prune 20 --blanklines 1 -N ")
#START SERVERS
print("Starting Servers...")
system(paste0(metamap.loc, "skrmedpostctl start"))
system(paste0(metamap.loc, "wsdserverctl start"))
#Wait 60 seconds before beginning run
Sys.sleep(60)

################### CONDUCT RUN ################### 
print("Beginning MetaMap Run...")
for(i in 1:length(input.files)) {
  if(i %% 2 == 0 & i != 1) { 
    print.noquote(paste0(i, " of ", length(input.files)))
    }
  system(paste0(metamap.base, input.files[i], " ", output.files[i]), 
        wait = T)
}
################### STOP SERVERS ################### 
print("Stopping Servers...")
system(paste0(metamap.loc, "skrmedpostctl stop"))
system(paste0(metamap.loc, "wsdserverctl stop"))



