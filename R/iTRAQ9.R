iTRAQ9 <- new("ReporterIons",
              description="8-plex iTRAQ and reporter + balance group",
              name="iTRAQ9",
              reporterNames=c(
                "iTRAQ9.113","iTRAQ9.114","iTRAQ9.115","iTRAQ9.116",
                "iTRAQ9.117","iTRAQ9.118","iTRAQ9.119","iTRAQ9.121",
                "iTRAQ9.305"),
              mz=c(113.1,114.1,115.1,116.1,
                117.1,118.1,119.1,121.1,
                305.1),
              col=c("yellow","blue","palegreen3","red",
                "tomato","black","purple","grey","orange"),
              width=0.05)
