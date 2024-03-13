
modelModules <- list()
modelModules$measurementModel <- 'f =~ X1 + X2 + X3 + X4 + X5 + X6'
modelModules$nuanceFactor <- 'n1 =~ X1 + X2 + X3

                              # uncorrelated factor and nuances
                              f ~~ 0*n1'
modelModules$allLatentY <- 'fY =~ 1 * y
                            y ~~ 0*y'


testSEMs <- list()

testSEMs$dgm <- paste0(modelModules$measurementModel, "\n",
                         modelModules$nuanceFactor, "\n",
                         modelModules$allLatentY, "\n",
                         '# structure model
                           fY ~ f + n1')

testSEMs$factor <- paste0(modelModules$measurementModel, "\n",
                           modelModules$allLatentY, "\n",
                           '# structure model
                             fY ~ f')

testSEMs$items <- paste0(modelModules$allLatentY, "\n",
                              '# structure model
                                fY ~  X1 + X2 + X3 + X4 + X5 + X6')

testSEMs$factorItems <- paste0(modelModules$measurementModel, "\n",
                                modelModules$allLatentY, "\n",
                                '# structure model
                                  fY ~ f + X1 + X2 + X3')

testSEMs$factorItem1 <- paste0(modelModules$measurementModel, "\n",
                                modelModules$allLatentY, "\n",
                                '# structure model
                                  fY ~ f + X1')


