iTRAQ9 <- new("ReporterIons",
              description="8-plex iTRAQ with isobaric tag",
              name="iTRAQ9",
              reporterNames=c(
                "iTRAQ9.113","iTRAQ9.114","iTRAQ9.115","iTRAQ9.116",
                "iTRAQ9.117","iTRAQ9.118","iTRAQ9.119","iTRAQ9.121",
                "iTRAQ9.305"),
              mz=c(113.13,114.13,115.13,116.13,
                117.13,118.13,119.13,121.13,
                305.13),
              col=c("yellow","blue","palegreen3","red",
                "tomato","black","purple","grey","orange"),
              width=0.05)
