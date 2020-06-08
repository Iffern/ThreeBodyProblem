#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(ggplot2)
library(gganimate)
library(gifski)
library(deSolve)


ui <- fluidPage(

    
    titlePanel("Three Body Problem"),

    
    sidebarLayout(
        sidebarPanel(
            
            h3("First body coefficients"),
            numericInput("m1",h4("Mass"),value=1),
            h4("Initial velocity"),
            splitLayout(numericInput("v1_x",h4("v1_x"),value=0,step=0.1),
                        numericInput("v1_y",h4("v1_y"),step=0.1,value=0)),
            h4("Initial position"),
            splitLayout(numericInput("x1",h4("x1"),value=-1,step=0.1),
                        numericInput("y1",h4("y1"),step=0.1,value=0)),
            h3("Second body coefficients"),
            numericInput("m2",h4("Mass"),value=1),
            h4("Initial velocity"),
            splitLayout(numericInput("v2_x",h4("v2_x"),value=0,step=0.1),
                        numericInput("v2_y",h4("v2_y"),step=0.1,value=0)),
            h4("Initial position"),
            splitLayout(numericInput("x2",h4("x2"),value=1,step=0.1),
                        numericInput("y2",h4("y2"),step=0.1,value=0)),
            h3("Third body coefficients"),
            numericInput("m3",h4("Mass"),value=1),
            h4("Initial velocity"),
            splitLayout(numericInput("v3_x",h4("v3_x"),value=0,step=0.1),
                        numericInput("v3_y",h4("v3_y"),step=0.1,value=0)),
            h4("Initial position"),
            splitLayout(numericInput("x3",h4("x3"),value=0,step=0.1),
                        numericInput("y3",h4("y3"),step=0.1,value=0)),
            numericInput("T",h3("Time"),value=120),
            numericInput("interval", h3("Time interval"),value=0.001, min=0, step = 0.001),
            numericInput("G",h3("Gravitational constant"), value=1,step=0.1),
            actionButton("start","Start")
        ),

        
        mainPanel(
           imageOutput("threeBodyPlot")
        )
    )
)


server <- function(input, output) {

    observeEvent(input$start, {
    output$threeBodyPlot <- renderImage({
        G <- input$G
        m1 <- input$m1
        m2 <- input$m2
        m3 <- input$m3
        params <- c(m1=m1,m2=m2,m3=m3,G=G)
        Y <- c(input$x1,input$y1,input$x2,input$y2,input$x3,input$y3,input$v1_x,input$v1_y,input$v2_x,input$v2_y,input$v3_x,input$v3_y)
        T <- seq(0,input$T,by=input$interval)
        ##H <- ((p1_x)^2+(p1_y)^2)/(2*m1) + ((p2_x)^2+(p2_y)^2)/(2*m2) + ((p3_x)^2+(p3_y)^2)/(2*m3) - G*m1*m2/sqrt((x1-x2)^2+(y1-y2)^2)- G*m1*m3/sqrt((x1-x3)^2+(y1-y3)^2)- G*m3*m2/sqrt((x3-x2)^2+(y3-y2)^2) 
        func <- function(t,y,params){
            with(as.list(c(params,y)), {
            dx1 <- y[7]
            dy1 <- y[8]
            dx2 <- y[9]
            dy2 <- y[10]
            dx3 <- y[11]
            dy3 <- y[12]
            dv1_x <- -G*m2*(y[1]-y[3])/((sqrt((y[1]-y[3])^2+(y[2]-y[4])^2))^3) - G*m3*(y[1]-y[5])/((sqrt((y[1]-y[5])^2+(y[2]-y[6])^2))^3)
            dv1_y <- -G*m2*(y[2]-y[4])/((sqrt((y[1]-y[3])^2+(y[2]-y[4])^2))^3) - G*m3*(y[2]-y[6])/((sqrt((y[1]-y[5])^2+(y[2]-y[6])^2))^3)
            dv2_x <- -G*m1*(y[3]-y[1])/((sqrt((y[1]-y[3])^2+(y[2]-y[4])^2))^3) - G*m3*(y[3]-y[5])/((sqrt((y[3]-y[5])^2+(y[4]-y[6])^2))^3)
            dv2_y <- -G*m1*(y[4]-y[2])/((sqrt((y[1]-y[3])^2+(y[2]-y[4])^2))^3) - G*m3*(y[4]-y[6])/((sqrt((y[3]-y[5])^2+(y[4]-y[6])^2))^3)
            dv3_x <- -G*m1*(y[5]-y[1])/((sqrt((y[1]-y[5])^2+(y[2]-y[6])^2))^3) - G*m2*(y[5]-y[3])/((sqrt((y[3]-y[5])^2+(y[4]-y[6])^2))^3)
            dv3_y <- -G*m1*(y[6]-y[2])/((sqrt((y[1]-y[5])^2+(y[2]-y[6])^2))^3) - G*m2*(y[6]-y[4])/((sqrt((y[3]-y[5])^2+(y[4]-y[6])^2))^3)
            list(c(dx1,dy1,dx2,dy2,dx3,dy3,dv1_x,dv1_y,dv2_x,dv2_y,dv3_x,dv3_y))
            })
        }
        df <- as.data.frame(ode(Y,T,func,params,method = "ode45"))
        outfile <- tempfile(fileext='.gif')
        anim = ggplot(data=df) +
                geom_point(mapping = aes(x= `1`,y=`2`,colour="red"),size=5) + 
                geom_point(mapping = aes(x= `3`,y= `4`,colour="blue"),size=5) + 
                geom_point(mapping = aes(x= `5`,y= `6`,colour="green"),size=5) + 
                geom_path(mapping = aes(x= `1`,y=`2`,colour = "red"))+
                geom_path(mapping = aes(x= `3`,y= `4`,colour="blue"))+
                geom_path(mapping = aes(x= `5`,y= `6`,colour="green"))+
                scale_colour_manual(values=c("red","blue","green"), labels=c("First body","Second body", "Third body"))+
                guides(color=guide_legend("Bodies"))+
                theme(legend.text = element_text(size=20), legend.title = element_text(size=30))+
                transition_reveal(along = time)
        
        anim_save("outfile.gif", animate(anim, height=700, width=900,renderer = gifski_renderer(loop=FALSE)))
        
        list(src = "outfile.gif",
             contentType = 'image/gif'
             # width = 400,
             # height = 300,
             # alt = "This is alternate text"
        )
    },deleteFile = TRUE)})
}

# Run the application 
shinyApp(ui = ui, server = server)
