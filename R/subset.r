## subset time frame data
Idents(obj.srt) <- 'time'
tmp <- subset(obj.srt, idents = c('naive','48hr_','72hr_') )
traj.set1 <- tmp
tmp <- subset(obj.srt, idents = c('72hr_','day5_','day30') )
traj.set2 <- tmp
tmp <- subset(obj.srt, idents = c('naive','day5_','day30') )
traj.set3 <- tmp

traj.set1 %>% saveRDS('2022.cell_hashing/imm_timecourse/traj.set1.rds')
traj.set2 %>% saveRDS('2022.cell_hashing/imm_timecourse/traj.set2.rds')
traj.set3 %>% saveRDS('2022.cell_hashing/imm_timecourse/traj.set3.rds')
