iTRAQ5 <- new("ReporterIons",
              description = "4-plex iTRAQ and reporter + balance group",
              name = "iTRAQ5",
              reporterNames = c("iTRAQ5.114", "iTRAQ5.115",
                                "iTRAQ5.116", "iTRAQ5.117",
                                "iTRAQ5.145"),
              mz = c(114.1112, 115.1083, 116.1116, 117.1150, 145.1),
              col = c("red", "green", "blue", "yellow", "grey"),
              width = 0.05)
