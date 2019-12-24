library(gdata)
library(ggpubr)

mut = read.xls('../results/cockroaches11_19/table.xlsx')

mut$Soldier = as.factor(mut$Soldier)

summary(mut$Worker)
summary(mut$Soldier)
summary(mut$Taxonomy)

workers = mut[mut$Worker == 1,]
workers = workers[!is.na(workers$Worker),]

a = as.data.frame(as.numeric(apply(workers[, 2:13], 2, sum)))
data = cbind(names(mut[, 2:13]), a)
names(data) = c('Subs', 'Number')

ggbarplot(data, 'Subs', 'Number', xlab="Substitution types", title = 'with workers', 
          fill = 'Subs', color = 'Subs', palette = c("#bdbdbd", "#7fcdbb", "#bdbdbd", "#bdbdbd", "#bdbdbd", "#feb24c", "#f03b20", "#bdbdbd", "#bdbdbd", "#bdbdbd", "#2c7fb8", "#bdbdbd"), combine = TRUE)

plot_hist = function(x, name){
  # x is dataframe with 12 column (subs types), name is a title of graph
  a = as.data.frame(as.numeric(apply(x, 2, sum)))
  data = cbind(names(x), a)
  names(data) = c('Subs', 'Number')
  
  ggbarplot(data, 'Subs', 'Number', xlab="Substitution types", title = name,
            fill = 'Subs', color = 'Subs', palette = c("#bdbdbd", "#7fcdbb", "#bdbdbd", "#bdbdbd", "#bdbdbd", "#feb24c", "#f03b20", "#bdbdbd", "#bdbdbd", "#bdbdbd", "#2c7fb8", "#bdbdbd"), combine = TRUE)
  
}

plot_hist(workers[, 2:13], 'workers')

withoutWorkers = mut[(mut$Worker == 0) & !is.na(mut$Worker),]

plot_hist(withoutWorkers[, 2:13], 'without workers')

summary(workers$Taxonomy)
table(withoutWorkers$Taxonomy)

####################################################################################
### Compare Ts/Tv ratio

# I want 2 density plots for each group

mut$TsTv = (mut$C_T + mut$A_G + mut$T_C + mut$G_A) / (mut$G_T + mut$A_T + mut$G_C +
                                                        mut$A_C + mut$C_A + mut$C_G +
                                                        mut$T_A + mut$T_G)

workers = mut[mut$Worker == 1,]
workers = workers[!is.na(workers$Worker),]
withoutWorkers = mut[(mut$Worker == 0) & !is.na(mut$Worker),]

filter_workers = rbind(workers, withoutWorkers)

summary(filter_workers$TsTv)

pdf('../results/cockroaches11_19/TsTsWorkersSoldiers.pdf')

ggplot(filter_workers, aes(TsTv, fill = filter_workers$Worker)) +
  # geom_histogram(aes(fill = filter_workers$Worker), alpha = 0.4) +
  geom_histogram( color="#e9ecef", alpha=0.4, position = 'dodge') +
  scale_fill_manual(values=c("#404080", "#69b3a2"))


# ggplot(filter_workers, aes(TsTv)) +
#   geom_density(aes(fill = filter_workers$Worker), alpha = 0.3)

# summary(mut$Soldier)
filter_soldiers = mut[!is.na(mut$Soldier),]

ggplot(filter_soldiers, aes(TsTv, fill = filter_soldiers$Soldier)) +
  # geom_histogram(aes(fill = filter_workers$Worker), alpha = 0.4) +
  geom_histogram( color="#e9ecef", alpha=0.4, position = 'dodge') +
  scale_fill_manual(values=c("#404080", "#69b3a2"))


dev.off()

summary(mut$Taxonomy)

cockroaches = c('Ectobiidae1', 'Tryonicidae', 'Blaberidae', 'Corydiidae', 'Ectobiidae2',
                'Lamproblattidae', 'Anaplectidae', 'Blattidae', 'Cryptocercidae', 'Ectobiidae3',
                'Nocticolidae')

mut = mut[!is.na(mut$Taxonomy),]

for(i in 1:nrow(mut)){
  if(mut$Taxonomy[i] %in% cockroaches){
    mut$Cockroaches[i] = 1
  }
  if(!(mut$Taxonomy[i] %in% cockroaches))
    {mut$Cockroaches[i] = 0}
}

mut$Cockroaches = as.factor(mut$Cockroaches)
summary(mut$Cockroaches)

pdf('../results/cockroaches11_19/TsTvCockroachesVsTermites.pdf')

ggplot(mut, aes(TsTv, fill = mut$Cockroaches)) +
  # geom_histogram(aes(fill = filter_workers$Worker), alpha = 0.4) +
  geom_histogram( color="#e9ecef", alpha=0.4, position = 'dodge') +
  scale_fill_manual(values=c("#69b3a2", "#404080"))

summary(mut$sumOfSubs)

plot_hist(mut[mut$Cockroaches == 1, 2:13], 'cockroaches')
plot_hist(mut[mut$Cockroaches == 0, 2:13], 'termites')

dev.off()
