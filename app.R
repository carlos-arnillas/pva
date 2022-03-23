#
# This is a Shiny web application. The core of the code can be extracted
# Most of the changes needed are highlighted below.
#

library(shiny)
library(data.table)
library(ggplot2)
library(gridExtra)

fun_txt0 <- "function(N, r0, K, sr, sk) {
    K <- K + rnorm(1,0,sk)
    if (K <= 50) K <- 50  # Never less than 50
    r <- r0*((K-N)/K) + rnorm(1,0,sr)
    Nt <- N + N*r
    if (Nt <= 0) {
        return(0)
    } else {
        return(Nt)
    }
}"


########
# The "ui" section can be deleted if not using shiny.
# Just remove the comments in the following lines and copy the section before this 
# and also the section below that is indicated as "the relevant code"
# ts <- 100          # time steps
# Ns <- numeric(ts)  # population
# reps <- 7
# No <- 10:10            # Initial population
# 
# 
# K <-  100          # Carrying capacity
# r0 <- 2.5
# s <- 0.1             # Growth stochasticity
# 
# Define UI for application that draws temporal trends in the population
ui <- fluidPage(

    # Application title
    titlePanel("Population viability analysis: Exploring a logistic model"),
    HTML("<p><i>Created by Carlos Alberto Arnillas</i></p>"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            numericInput("ts",
                        "Years [t] (1 to 1000):",
                        min = 1,
                        max = 1000,
                        value = 100),
            numericInput("reps",
                         "Repetitions (1 to 500):",
                         min = 1,
                         max = 500,
                         value = 1),
            h5("Initial population [No] (1 to 1000)"),
            splitLayout(numericInput("No1",
                         "From:",
                         min = 1,
                         max = 1000,
                         value = 50),
                    numericInput("No2",
                         "To:",
                         min = 1,
                         max = 1000,
                         value = 50)),
            numericInput("sr",
                         "Growth-related noise [s] (0 to 5):",
                         min = 0,
                         max = 1,
                         value = 0),
            numericInput("K",
                         "Carrying capacity [K] (50 to 1000):",
                         min = 1,
                         max = 1000,
                         value = 500),
            numericInput("sk",
                         "Carrying capacity-related noise [sk] (0 to 1000):",
                         min = 0,
                         max = 1000,
                         value = 0),
            numericInput("r0",
                         "Maximum population growth rate at low density [r0] (-5 to 5):",
                         min = -5,
                         max = 5,
                         value = 0.5),
            textInput("lcuts",
                      "Check points (separated by commas):",
                      value = "10,30,50,100"),
            textAreaInput("fun_txt",
                          "Function (disabled):",
                          height="250px",
                          value = fun_txt0)
        ),

        # Show a plot of the generated distribution
        mainPanel(
            tabsetPanel(
                tabPanel(
                    "Plots",
                    h3("Population trend"),
                    HTML("<p>Each line represents a repetition (potential history), coloured according its initial population. The black line is one repetition, coloured in black to help visualization, while the red dotted line indicates the carring capacity (<i>K</i>).</p>"),
                    HTML("<p>Check points are represented with vertical dashed lines and the corresponding number indicates how many populations are alive (<i>N<sub>t</sub></i>&nbsp;>&nbsp;0) that year.</p>"),
                    plotOutput("distPlot1", height="200px"),
                    h3("Final conditions"),
                    p("Number of population histories (lines) with a certain abundance at the end of the simulation. The red dotted line indicates the carring capacity (K)."),
                    plotOutput("distPlot2", height="200px")
                ),
                tabPanel(
                    "Table",
                    p("Values correspond to the series used to plot the black line."),
                    tableOutput("tblTrends")
                )
            )
        )
    )
)

validNum <- function(x, min, max) {
    if (is.na(x)) return(FALSE)
    return(!(is.na(as.numeric(x)) | x < min | x > max))
}

# Define server logic required to draw temporal trends
server <- function(input, output) {
    
    resF <- reactive({
        ####### All the relevant code starts here!
        # Interpreting the function: fixed for now
        # r <- eval(parse(text=input$fun_txt))
        error <- ""
        r <- eval(parse(text=fun_txt0))
        # loading the variables
        ts <- input$ts
        reps <- input$reps
        K <- input$K
        sk <- input$sk
        r0 <- input$r0
        sr <- input$sr
        No1 <- input$No1
        No2 <- input$No2
        # print("xxx")
        # print(str(as.numeric(ts)))
        # print("x-x")
        # print(is.na(ts))
        # Checking values
        #pass <<- FALSE
        if (!validNum(ts, 10, 1000)) error <- "Adjust number of years!"
        if (!validNum(reps, 1, 500)) error <- "Adjust repetitions!"
        if (!validNum(No1, 1, 1000)) error <- "Adjust initial population (from)!"
        if (!validNum(No2, 1, 1000)) error <- "Adjust initial population (to)!"
        if (!validNum(K, 50, 1000)) error <- "Adjust carrying capacity!"
        if (!validNum(sk, 0, 1000)) error <- "Adjust carrying capacity variability!"
        if (!validNum(r0, -5, 5))  error <- "Adjust carrying intrinsic growth rate!"
        if (!validNum(sr, 0, 5)) error <- "Adjust intrinsic growth rate variability!"

        if (error != "") {
            error 
        } else {
            # Let's proceed
            Ns <- numeric(ts)  # population
            No <- input$No1:input$No2
            # Accommodate No and reps (reps drives the number of cases to explore)
            if (length(No) != reps) {
                Nox <- numeric(reps)
                Nox[] <- No
                No <- Nox
                rm(Nox)
            }
            
            # our_seed <- 78L
            # set.seed(our_seed)
            
            # Prep the output table
            res <- data.table(reps.=integer(0), No.=numeric(0), ts.=integer(0), Ns.=numeric(0))
            # Running the model
            for (k in 1:reps) {
                Ns[1] <- No[k]
                for (i in 1:(ts-1)) {
                    Ns[i+1] <- r(Ns[i], r0, K, sr, sk)
                }
                # Storing the data
                res <- rbind(res, cbind(reps.=k, No.=No[k], ts.=1:ts, Ns.=Ns))
            }
            res
        }
    })
    
    # Checkpoints
    lcutsF <- reactive({        
        lcuts <- input$lcuts
        if (lcuts=="") {
            lcuts <- NA
        } else if (!grepl("^[0-9]+(,[0-9]+)*$", lcuts)) error <- "Adjust check points!"
        if (!is.na(lcuts)) {
            lcuts <- sort(unique(as.integer(strsplit(lcuts, ",", fixed=TRUE)[[1]])))
            lcuts <- lcuts[lcuts > 0]
            lcuts <- lcuts[lcuts <= input$ts]
        }
        lcuts
    })

    output$distPlot1 <- renderPlot({
        res <- resF()
        # Any errors?
        if (is.character(res)) {
            ggplot() + annotate("text", x = 4, y = 25, size=8, label = res) + theme_void()
        } else {
            # Plotting
            survived <- nrow(resF()[Ns. > 0 & ts.==max(ts.)])
            q1 <- ggplot(resF(), aes(x=ts., y=Ns.)) + 
                geom_hline(yintercept=input$K, colour="red", linetype="dotted") +
                geom_line(aes(group=reps., colour=factor(No.)), alpha=0.5) + 
                geom_line(data=resF()[reps.==min(reps.)]) + 
                labs(x=sprintf("Year\n\nSurviving populations: %d/%d", survived, input$reps), colour="No") + theme_bw()
            lcutsX <- lcutsF()
            if (is.character(lcutsX)) {
                q1 <- q1 + annotate(x=0, y=0, label=lcutsX, colour="gray")
            } else {
                if (!is.na(lcutsX[1])) {
                    cuts <- resF()[Ns. > 0 & ts. %in% lcutsX, by=ts., .(n.=.N)]
                    q1 <- q1 + geom_vline(data=cuts, aes(xintercept=ts.), linetype="dashed", colour="gray70") +
                          geom_label(data=cuts, y=Inf, vjust="inward", hjust=0.5, aes(label=n.))
                }
            }
            # grid.arrange(q1,q2,ncol=1,heights=list(90, 50))
            q1
    
            ####### All the relevant code finishes here!
        }
    })
    output$distPlot2 <- renderPlot({
        res <- resF()
        # Any errors?
        if (is.character(res)) {
            ggplot() + annotate("text", x = 4, y = 25, size=8, label = res) + theme_void()
        } else {
            
            q2 <- ggplot(resF()[ts.==max(ts.)], aes(x=Ns.)) + geom_histogram(bins=50) + 
                coord_cartesian(xlim=c(if (length(unique(resF()[ts.==max(ts.),Ns.]))==1) NA else 0, NA)) + 
                labs(x="Final population") + 
                geom_vline(xintercept=input$K, colour="red", linetype="dotted", alpha=0.5) + theme_bw()
            q2
        }
    })
    
    output$tblTrends <- renderTable({
        res <- resF()
        if (is.character(res)) {
            res <- data.table(`Error:`=res)
        } else {
            res <- res[reps.==min(reps.)][,.(year=ts., N=Ns.)]
        }
        res
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
