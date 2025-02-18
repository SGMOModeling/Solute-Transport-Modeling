library(data.table)
library(plotly)
library(reshape)
library(contourPlot)
library("colorspace")
library(grDevices)
library(RColorBrewer)
#library(tidyverse)
library(akima)
library(magick)


library(sf)
library(raster)

# Number of contour intervals
NLs=5  

NL_flood=20



dx=300
dy=300

iskip=10   # output skip in Fortran


#var = readline(prompt = "Enter the layer number : ")

#ilyr= as.integer(var);

#Lower left coordinates
# XLL=620000.0
# YLL=4247150.0 
# 
# XUR=660000.0
# YUR=4280043.0


XLL=660000
YLL= 4201234 

XUR=675000
YUR= 4220000

ilyr=1     # XY conc
irow=31    # xZ conc
icol=36   # YZ conc



#path="C:/Users/ubandara/Documents/TransportModels/AdvectionDiffusion/Groundwater/3-DCase/AdvDisp3DLagrPar/Output/C2VSimFG/HypotheticalCase/Output/ReleaseLyr2/Bilinear_Linear_Comb/"
#path="C:/Users/ubandara/Documents/SoluteTransport/Models/3D/AdvDisp3DLagrPar/Output/C2VSimFG/HypotheticalCase/Output/ReleaseLyr2/Bilinear_Linear_Comb/"

path="C:/Users/ubandara/Documents/SoluteTransport/Models/3D/AdvDisp3DLagrPar/Output/"

# Load shapefile (replace with your own shapefile path)
shapefile_path <- paste0(path,"C2VSimFG_Elements.shp")
shape_data <- st_read(shapefile_path)

#Prepare the same color ramp
lvls <-seq(from = 0, to = 10, by =0.1)

file_name=paste("ConLagrng_3D_AllInOne.dat", sep="")    
ConcLag_table<-fread(paste0(path,file_name), header="auto")

ConcLag_table<-t(ConcLag_table)

#Convert the table into data frame
ConcLag<-as.data.frame(ConcLag_table, stringsAsFactors = FALSE)

#Use the first row as the column titles
colnames(ConcLag)<-ConcLag[1,] # First row becomes column names
ConcLag<-ConcLag[-1,]  # Now remove the first row

ConcLag[] <- lapply(ConcLag, function(x) as.numeric(as.character(x)))

#Extact times column headers dynamically starts with "Conc_"

TimeSeries<-grep("^Conc_", colnames(ConcLag), value=TRUE) 


  # Loop through TimeSeries to access columns in XY_ConcLag


NumFigures=42

#for (i in seq_along(TimeSeries)) {
istart_fig=3
  for (i in istart_fig:NumFigures) {
    # Construct the column name
    col_name <-TimeSeries[i]
    
    # Access the column using the constructed name
    if (col_name %in% colnames(ConcLag)) { # Ensure the column exists
      conc_data <-(ConcLag[[col_name]])
      print(paste("Processing column:", col_name))
  
          #XY Plot
          fig_name = paste("Lagrangian_ColorFlood_XY_Lyr_",ilyr,"_", i, ".png", sep="")
          figure<-png(paste0(path,fig_name),type="cairo")
          
          #select specified layer number to plot
          XY_ConcLag <- ConcLag[ConcLag$LAYER==ilyr, ]
          
          # Remove duplicates from data
          XY_ConcLag <- XY_ConcLag %>%
            distinct(XCOORD, YCOORD, .keep_all = TRUE)
          
          # Extract the necessary columns
          x <- XY_ConcLag$XCOORD
          y <- XY_ConcLag$YCOORD
          XYConc <-XY_ConcLag[[TimeSeries[i]]] 
          
          # Remove data points with zero values
          non_zero_indices <- XYConc != 0  # Logical vector identifying non-zero values
          x_non_zero <- x[non_zero_indices]
          y_non_zero <- y[non_zero_indices]
          XYConc_non_zero <- XYConc[non_zero_indices]
          
          
          # Interpolate data
          #my.matrix <- interp(x, y, XYConc)
          my.matrix <- interp(x_non_zero, y_non_zero, XYConc_non_zero)
          
          # Check for and handle any NA values in the interpolated results
          my.matrix$z[is.na(my.matrix$z)] <- 0  # Replace NAs with 0 or another value if necessary
          
          # Define levels for contouring, you can customize this
          #lvls <- pretty(my.matrix$z, n = 20)
          
          
          # Create the filled contour plot
          # filled.contour(my.matrix,
          #                color.palette = colorRampPalette(c('white','blue','green','yellow','red','darkred')),
          #                levels = lvls,
          #                asp = 1,
          #                xlim = c(622000.0, 626000.0), 
          #                ylim = c(4257000,4261000),
          #                plot.title = {par(cex.main = 1); title(main = paste("Lagrangian Particle Solution for 3D Hypothetical Case-C2VSimFG \n in XY plane Layer", ilyr),
          #                                                       xlab = "Distance along X-axis", ylab = "Distance along Y-axis")},
          #                key.title = {par(cex.main = 0.7); title(main = "Concentration\n(mg/l)")})
          # mtext(paste("Month:", i*iskip), side = 1, line = -23, adj = 0.6, cex = .9)
          # dev.off()
          
          # # Add Mesh
          #Adjust the margins
          par(mar = c(4, 4, 4, 4))  # Adjust values as needed: bottom, left, top, right

          # Create the filled contour plot
          filled.contour(
            my.matrix,
            color.palette = colorRampPalette(c('white','blue','green','yellow','red','darkred')),
            levels = lvls,
            asp = 1,
            xlim = c(XLL, XUR),
            ylim = c(YLL, YUR),

            # Title and axis labels
            plot.title = {
              #par(cex.main = 1, mgp = c(5, 1, 0), las = 3)
              par(cex.main = 1)  # Change `las` as desired
              title(main = paste("Lagrangian Parcel Solution for 3D Hypothetical Case \n C2VSimFG: XY plane Layer", ilyr),
                    xlab = "X(m)", ylab = "Y(m)")
            },

            # Key title
            key.title = {
              par(cex.main = 0.7)
              title(main = "Concentration\n(mg/l)")
            },

            # Overlay shapefile on the contour plot
            plot.axes = {
              par(mgp = c(2.5, 1, 0), las = 0)  # Increase y-axis title distance
              axis(1); axis(2)
              plot(st_geometry(shape_data), add = TRUE, col = "transparent", border = "grey", lwd = 1)
            }
          )
          mtext(paste("Month:", i*iskip), side = 1, line = -23, adj = 0.6, cex = .9)
          dev.off()
          # 
          
          
          
    }
          
       print(paste("Processing XZ Plot:"))
      # # XZ plot
      fig_name = paste("Lagrangian_ColorFlood_XZ_Col_",icol,"_", i, ".png", sep="")
      figure<-png(paste0(path,fig_name),type="cairo")

      #select specified column number to plot
      #XZ_ConcLag <- ConcLag[ConcLag$YCOORD==irow*dy+YLL,  ]
      XZ_ConcLag <- ConcLag[ConcLag$COL==icol, ]

      # Remove duplicates from the X and Z coordinates
      XZ_ConcLag <- XZ_ConcLag %>%
        distinct(XCOORD, ZCOORD, .keep_all = TRUE)

      # Extract the necessary columns for XZ view
      x <- XZ_ConcLag$XCOORD  # X-coordinate
      z <- XZ_ConcLag$ZCOORD  # Z-coordinate
      XZConc <- XZ_ConcLag[[TimeSeries[i]]]   # Concentration

      # Remove data points with zero values
      non_zero_indices <- XZConc != 0  # Logical vector identifying non-zero values
      x_non_zero <- x[non_zero_indices]
      z_non_zero <- z[non_zero_indices]
      XZConc_non_zero <- XZConc[non_zero_indices]


      # Interpolate data
      my.matrix <- interp(x_non_zero, z_non_zero, XZConc_non_zero)


      # Check for and handle any NA values in the interpolated results
      my.matrix$z[is.na(my.matrix$z)] <- 0  # Replace NAs with 0 or another value if necessary

      # Define levels for contouring, customize this as needed
      #lvls <- pretty(my.matrix$z, n = 20)


      filled.contour(my.matrix, color.palette=colorRampPalette(c('white','blue','green','yellow','red','darkred')),levels=lvls,asp=1.0,
                     xlim = range(660000, 670000.0),
                     ylim = range(-1500:500),
                     plot.title = {par(cex.main=1);title(main = paste("Lagrangian Particle Solution for 3D Hypothetical Case-C2VSimFG \n in XZ plane Column",icol),
                                                         xlab = "Distance along X-axis", ylab = "Distance along Z-axis")},
                     key.title = {par(cex.main=0.7);title(main="Concentration\n(mg/l)")})
      mtext(paste("Month:", i*iskip), side = 1, line = -23, adj = 0.6, cex = .9)
      dev.off()

      #

      print(paste("Processing YZ Plot:"))
      # YZ plot
      fig_name = paste("Lagrangian_ColorFlood_YZ_Row_",irow,"_", i, ".png", sep="")
      figure<-png(paste0(path,fig_name),type="cairo")

      #select specified column number to plot
      #YZ_ConcLag <- ConcLag[ConcLag$XCOORD==icol*dx+XLL,  ]
      YZ_ConcLag <- ConcLag[ConcLag$ROW==irow, ]


      # Remove duplicates from data
      YZ_ConcLag <- YZ_ConcLag %>%
        distinct(YCOORD, ZCOORD, .keep_all = TRUE)

      # Extract the necessary columns
      y <- YZ_ConcLag$YCOORD
      z <- YZ_ConcLag$ZCOORD
      YZConc <- YZ_ConcLag[[TimeSeries[i]]]


      # Remove data points with zero values
      non_zero_indices <- YZConc != 0  # Logical vector identifying non-zero values
      y_non_zero <- y[non_zero_indices]
      z_non_zero <- z[non_zero_indices]
      YZConc_non_zero <- YZConc[non_zero_indices]


      # Interpolate data
      my.matrix <- interp(y_non_zero, z_non_zero, YZConc_non_zero)

      # Interpolate data
      #my.matrix <- interp(x, y, YZConc)

      # Check for and handle any NA values in the interpolated results
      my.matrix$z[is.na(my.matrix$z)] <- 0  # Replace NAs with 0 or another value if necessary

      # Define levels for contouring, you can customize this
      #lvls <- pretty(my.matrix$z, n = 20)



      filled.contour(my.matrix, color.palette=colorRampPalette(c('white','blue','green','yellow','red','darkred')),levels=lvls,
                     asp=1.0,xlim = range(4200000.0, 4215000.0),
                     ylim = range(-1500:500),
                     plot.title = {par(cex.main=1);title(main = paste("Lagrangian Particle Solution for 3D Hypothetical Case-C2VSimFG \n in YZ plane row",irow),
                                                         xlab = "Distance along Y-axis", ylab = "Distance along Z-axis")},
                     key.title = {par(cex.main=0.7);title(main="Concentration\n(mg/l)")})
      mtext(paste("Month:", i*iskip), side = 1, line = -23, adj = 0.6, cex = .9)
      dev.off()

}



# #Create an animation  (use magick library)
frames <- paste0(path, "Lagrangian_ColorFlood_XY_Lyr_",ilyr,"_",istart_fig:NumFigures, ".png")
m <- image_read(frames)
m <- image_animate(m, fps = 4, loop = 0)
#image_write(m, "RMA_Contaminant_Transport_Lagrn.gif")
image_write(m, paste0(path,"XY_view_of_C2VSimFG_Hypothetical_Transport_Lagrang_Lyr_",ilyr,"_", NumFigures*iskip, "Months.gif"))


frames <- paste0(path, "Lagrangian_ColorFlood_XZ_Col_",icol,"_",istart_fig:NumFigures, ".png")
m <- image_read(frames)
m <- image_animate(m, fps = 4, loop = 0)
#image_write(m, "RMA_Contaminant_Transport_Lagrn.gif")
image_write(m, paste0(path,"XZ_view_of_C2VSimFG_Hypothetical_Transport_Lagrang_Row_",irow,"_", NumFigures*iskip, "Months.gif"))


#
frames <- paste0(path, "Lagrangian_ColorFlood_YZ_Row_",irow,"_",istart_fig:NumFigures, ".png")
m <- image_read(frames)
m <- image_animate(m, fps = 4, loop = 0)
#image_write(m, "RMA_Contaminant_Transport_Lagrn.gif")
image_write(m, paste0(path,"YZ_view_of_C2VSimFG_Hypothetical_Transport_Lagrang_Col_",icol,"_", NumFigures*iskip, "Months.gif"))
#


















#Final color flood plots in the same figure------------------------------------------------------------------------------------------!




#Prepare the same color ramp
lvls <-seq(from =1, to = 11, by =0.25)

fontsize=1.2

  col_name <-TimeSeries[i]

    conc_data <-(ConcLag[[col_name]])
    print(paste("Processing column:", col_name))
    
    #XY Plot
    fig_name = paste("Lagrangian_ColorFlood_35_Yr_XY_Lyr_",ilyr,"_", i, ".png", sep="")
    figure<-png(paste0(path,fig_name),type="cairo")
    

    #select specified layer number to plot
    XY_ConcLag <- ConcLag[ConcLag$LAYER==ilyr, ]
    
    # Remove duplicates from data
    XY_ConcLag <- XY_ConcLag %>%
      distinct(XCOORD, YCOORD, .keep_all = TRUE)
    
    # Extract the necessary columns
    x <- XY_ConcLag$XCOORD
    y <- XY_ConcLag$YCOORD
    XYConc <-XY_ConcLag[[TimeSeries[i]]] 
    
    # Remove data points with zero values
    non_zero_indices <- XYConc != 0  # Logical vector identifying non-zero values
    x_non_zero <- x[non_zero_indices]
    y_non_zero <- y[non_zero_indices]
    XYConc_non_zero <- XYConc[non_zero_indices]
    
    
    # Interpolate data
    #my.matrix <- interp(x, y, XYConc)
    my.matrix <- interp(x_non_zero, y_non_zero, XYConc_non_zero)
    
    # Check for and handle any NA values in the interpolated results
    my.matrix$z[is.na(my.matrix$z)] <- 0  # Replace NAs with 0 or another value if necessary

    
    # Create the filled contour plot
      #Adjust the margins
    par(mar = c(4, 4, 4, 4))  # Adjust values as needed: bottom, left, top, right
    
    filled.contour(
      my.matrix,
      color.palette = colorRampPalette(c('white','blue','green','yellow','red','darkred')),
      levels = lvls,
      asp = 1,
      xlim = c(663000.0, 668000.0), 
      ylim = c(4206000.0, 4214000.0),
      
      # Title and axis labels
      plot.title = {
        par(cex.main=1,cex.lab = 1.2); 
        title(main = paste(""),
        xlab = "X(m) ", ylab = "Y(m)")
        },
      key.title = {
        par(cex.main = 0.8);
        title(main = "Conc.\n(mg/l)")
        },
      # Overlay shapefile on the contour plot
      plot.axes = {
        par(mgp = c(2.5, 1, 0), las = 0)  # Increase y-axis title distance
        axis(1,las = 0, cex.axis = 1.2); axis(2,las = 0, cex.axis = 1.2)
        plot(st_geometry(shape_data), add = TRUE, col = "transparent", border = "grey", lwd = 1)
        #grid()
        #grid(nx = 100, ny = 100, col = "lightgray", lty = "dotted", lwd = 0.5) 
        points(665786,  4210606, pch=18, col="magenta", bg="blue", cex=2.0)
        # Add the symbol for the legend
        points(664800,4206500, pch=18, col="magenta", bg="yellow", cex=2.0)
      }
    )
    
    
    mtext(paste("(a) XY view through layer ", ilyr), side = 1, line = -23, adj = 0.05, cex = fontsize)

    mtext(paste("Release Location"), side = 1, line = -2.0, adj = 0.01, cex = 1.0)
    
    dev.off() 
    
    
  
    asp_ratio=2  
  
  # # XZ plot
    fig_name = paste("Lagrangian_ColorFlood_35_yr_XZ_Row_",icol,"_", i, ".png", sep="")
    figure<-png(paste0(path,fig_name),type="cairo")


  #select specified column number to plot
  XZ_ConcLag <- ConcLag[ConcLag$COL==icol, ]

  # Remove duplicates from the X and Z coordinates
  XZ_ConcLag <- XZ_ConcLag %>%
    distinct(XCOORD, ZCOORD, .keep_all = TRUE)

  # Extract the necessary columns for XZ view
  x <- XZ_ConcLag$XCOORD  # X-coordinate
  z <- XZ_ConcLag$ZCOORD  # Z-coordinate
  XZConc <- XZ_ConcLag[[TimeSeries[i]]]   # Concentration

  # Remove data points with zero values
  non_zero_indices <- XZConc != 0  # Logical vector identifying non-zero values
  x_non_zero <- x[non_zero_indices]
  z_non_zero <- z[non_zero_indices]
  XZConc_non_zero <- XZConc[non_zero_indices]


  # Interpolate data
  my.matrix <- interp(x_non_zero, z_non_zero, XZConc_non_zero)


  # Check for and handle any NA values in the interpolated results
  my.matrix$z[is.na(my.matrix$z)] <- 0  # Replace NAs with 0 or another value if necessary


  filled.contour(my.matrix, color.palette=colorRampPalette(c('white','blue','green','yellow','red','darkred')),levels=lvls,asp=asp_ratio,
                 xlim = c(664000.0, 667500.0),
                 ylim = range(-1400:0),
                 plot.title = {par(cex.main=1,cex.lab = 1.2);title(main = paste(""),
                                                     xlab = "X(m)", ylab = "Z(m)")},
                 key.title = {par(cex.main=0.8);title(main="Conc.\n(mg/l)")},
                 plot.axes = {
                   # Add x-axis with labels parallel to the axis
                   axis(1, las = 0, cex.axis = 1.2)  # las = 0 makes labels parallel to x-axis
                   # Add y-axis with labels parallel to the axis
                   axis(2, las = 0, cex.axis = 1.2)  # las = 0 makes labels parallel to y-axis
                   cex.axis = 1.5# Optionally, add gridlines or points
                   points(665786, -100, pch=18, col="magenta", bg="blue", cex=2.0)
                   # Add the symbol for the legend
                   points(665400,-1570, pch=18, col="magenta", cex=2.0)
                   grid()
                 },)


  mtext(paste("(b) XZ view through column :",icol), side = 1, line = -23, adj = 0.05, cex = fontsize)
  mtext(paste("Vertical Exaggeration:",asp_ratio ), side = 1, line = -2.0, adj = 0.55, cex = 1.0)
  mtext(paste("Release Location"), side = 1, line = -2.0, adj = 0.01, cex = 1.0)
  #mtext(paste("★"), side = 1, line = -2.0, adj = 0.3, cex = 1.0)
 dev.off()




  # YZ plot

 fig_name = paste("Lagrangian_ColorFlood_35_Yr_YZ_Col_",irow,"_", i, ".png", sep="")
 figure<-png(paste0(path,fig_name),type="cairo")

  #select specified column number to plot
  YZ_ConcLag <- ConcLag[ConcLag$ROW==irow, ]


  # Remove duplicates from data
  YZ_ConcLag <- YZ_ConcLag %>%
    distinct(YCOORD, ZCOORD, .keep_all = TRUE)

  # Extract the necessary columns
  y <- YZ_ConcLag$YCOORD
  z <- YZ_ConcLag$ZCOORD
  YZConc <- YZ_ConcLag[[TimeSeries[i]]]


  # Remove data points with zero values
  non_zero_indices <- YZConc != 0  # Logical vector identifying non-zero values
  y_non_zero <- y[non_zero_indices]
  z_non_zero <- z[non_zero_indices]
  YZConc_non_zero <- YZConc[non_zero_indices]


  # Interpolate data
  my.matrix <- interp(y_non_zero, z_non_zero, YZConc_non_zero)


  # Check for and handle any NA values in the interpolated results
  my.matrix$z[is.na(my.matrix$z)] <- 0  # Replace NAs with 0 or another value if necessary


  filled.contour(my.matrix, color.palette=colorRampPalette(c('white','blue','green','yellow','red','darkred')),levels=lvls,
                 asp=asp_ratio,xlim = c(4208000.0, 4212000.0),
                 ylim = range(-1300:0),
                 plot.title = {par(cex.main=1,cex.lab = 1.2);title(main = paste(""),
                                                     xlab = "Y(m)", ylab = "Z(m)")},
                 key.title = {par(cex.main=1.0);title(main="Conc.\n(mg/L)")},
                 plot.axes = {
                   # Add x-axis with labels parallel to the axis
                   axis(1, las = 0, cex.axis = 1.2)  # las = 0 makes labels parallel to x-axis
                   # Add y-axis with labels parallel to the axis
                   axis(2, las = 0, cex.axis = 1.2)  # las = 0 makes labels parallel to y-axis
                   # Optionally, add gridlines or points
                   points(4210606, -100, pch=18, col="magenta", bg="blue", cex=2.0)
                   # Add the symbol for the legend
                   points(4209600,-1650, pch=18, col="magenta", cex=2.0)
                   grid()
                 },)


  mtext(paste("(c) YZ view through row:", irow), side = 1, line = -23, adj = 0.05, cex = fontsize)
  mtext(paste("Vertical Exaggeration:",asp_ratio ), side = 1, line = -2.0, adj = 0.55, cex = 1.0)
  mtext(paste("Release Location"), side = 1, line = -2.0, adj = 0.01, cex = 1.0)


dev.off() 


#End of Final plot in the same figure------------------------------------------------------------------------------------------!


















#--Create Contour Plots---------------------------------------------------------------------------------------------------------------------------------------------------!
#XY
fig<-tiff(paste0(path,"Lagrng_Contour_XY_Lyr_",ilyr,".tiff"),type="cairo", res=100)


lvls <-c( 5, 10,20,30)
# Create the contour plot
x <- XY_ConcLag$XCOORD
y <- XY_ConcLag$YCOORD
XYConc <- XY_ConcLag$CONC_PPM
my.matrix  <- interp(x,y,XYConc)
contour(my.matrix, levels=lvls, asp=1,
        xlim = range(622000.0, 626000.0, finite = TRUE),
        ylim = range(4257000,4261000, finite = TRUE),
        col='green',
        vfont = c("sans serif", "bold"), 
        labcex = 1.0,
        #cex.lab = 2.5,cex.axis = 3.0,
        lwd = 2, 
        lty = 1, 
        xlab = "X(m)", 
        ylab = "Y(m)")

# Add single star at specific coordinate
points(563019.0, 4429475, pch=8, col="red", cex=1.0)

mtext(paste("(a) XY View in Layer",ilyr), side = 1, line = -14, adj = 0.1, cex = 0.8)
mtext(paste("Vertical Exaggeration :",1), side = 1, line = -1.5, adj = 0.9, cex = 0.8)
dev.off() 

#XZ
fig<-tiff(paste0(path,"Lagrng_Contour_XZ_Row_",irow,".tiff"),type="cairo", res=100)


#lvls <-c( 10,20,30)
asp_ratio=0.1
x<-XZ_ConcLag$XCOORD
z<-XZ_ConcLag$ZCOORD
my.matrix  <- interp(x,z,XZConc)
contour(my.matrix, levels=lvls,xlim = range(622000.0, 626000.0, finite = TRUE),ylim = range(-800,500, finite = TRUE),col='green',asp=asp_ratio,labels =lvls,
        vfont = c("sans serif", "bold"), labcex = 1.0,cex.lab = 1.0,cex.axis = 1.0,lwd = 2, lty = 1, xlab = "X(m)", ylab = "Z(m)")


#legend(x='topright',legend="",
#       col=c("green","red"), pch=c(NA,NA),lty = c(2),  cex=0.8,lwd = 2)
#mtext(paste("Days= 100"), side = 1, line = -16, adj = 0.1, cex = 1.0)
#mtext(paste("Pe=1"), side = 1, line = -14, adj = 0.1, cex = 0.8)   
mtext(paste("(b) XZ view through row",irow), side = 1, line = -14, adj = 0.1, cex = 0.8)
mtext(paste("Vertical Exaggeration :",asp_ratio), side = 1, line = -1.5, adj = 0.9, cex = 0.8)
dev.off() 

#YZ
fig<-tiff(paste0(path,"Lagrng_Contour_YZ_Col_",icol,".tiff"),type="cairo", res=100)

#lvls <-c( 10,20,30)

asp_ratio=0.1
y<-YZ_ConcLag$YCOORD
z<-YZ_ConcLag$ZCOORD
my.matrix  <- interp(y,z,YZConc)
contour(my.matrix, levels=lvls,xlim = range(4257000,4261000, finite = TRUE),ylim = range(-800,500, finite = TRUE),col='green',asp=asp_ratio,
        vfont = c("sans serif", "bold"), labcex = 1.0,cex.lab =1.0,lwd = 2,cex.axis = 1.0, lty = 1, xlab = "Y(m)", ylab = "Z(m)")

points( 4429475,-50, pch=8, col="red", cex=1.0)
#legend(x='topright',legend="",
#       col=c("green","red"), pch=c(NA,NA),lty = c(2),  cex=0.8,lwd = 2)
#mtext(paste("Days= 100"), side = 1, line = -16, adj = 0.1, cex = 1.0)
#mtext(paste("Pe=1"), side = 1, line = -14, adj = 0.1, cex = 0.8)   
mtext(paste("(c) YZ view through column",icol), side = 1, line = -14, adj = 0.1, cex = 0.8)
mtext(paste("Vertical Exaggeration:",asp_ratio ), side = 1, line = -1.5, adj = 0.9, cex = 0.8)
dev.off() 









# Plot all three in one figure-------------------------------------------------------------------------------#
# Set up the plotting layout: 1 row, 3 columns

#fig<-tiff(paste0(path,"Lagrng_Contour_",ilyr,"_",irow,"_",icol,".tiff"),type="cairo", res=100)

fig<-png(paste0(path, "Lagrng_Contour_", ilyr, "_", irow, "_", icol, ".png"), width = 1800, height = 600)

par(mfrow = c(1, 3))  

# Ensure sufficient margins (optional)
par(mar = c(8, 8, 4,3))  # Adjust margins to fit multiple plots

# Increase mgp[2] to 1.5 for more space between axis and tick labels
par(mgp = c(6, 3.0, 0)) 

#XY
lvls <-c(1,5, 10,20,30)
# Create the contour plot
x <- XY_ConcLag$XCOORD
y <- XY_ConcLag$YCOORD
XYConc <- XY_ConcLag$CONC_PPM
my.matrix  <- interp(x,y,XYConc)
contour(my.matrix, levels=lvls, asp=1,
        xlim = range(560000:580000, finite = TRUE),
        ylim = range(4420000:4440000, finite = TRUE),
        col='green',
        vfont = c("sans serif", "bold"), 
        labcex = 2.0,
        cex.lab = 2.5,cex.axis = 3.0,
        lwd = 2, 
        lty = 1, 
        xlab = "X(m)", 
        ylab = "Y(m)")

# Add single star at release location
points(563019.0, 4429475, pch=17, col="red", cex=3.0)

# Add the symbol for the legend
points(566500.0, 4420000, pch=17, col="red", cex=3.0)

mtext(paste("(a) XY View in Layer",ilyr), side = 1, line = -45, adj = 0.1, cex = 2.0)
mtext(paste("Vertical Exaggeration :",1), side = 1, line = -2.0, adj = 0.9, cex = 1.5)
mtext(paste("Release location"), side = 1, line = -2.0, adj = 0.05, cex = 1.5)
#dev.off() 

#XZ
#lvls <-c( 10,20,30)
asp_ratio=10
x<-XZ_ConcLag$XCOORD
z<-XZ_ConcLag$ZCOORD
my.matrix  <- interp(x,z,XZConc)
contour(my.matrix, levels=lvls,xlim = range(560000,580840, finite = TRUE),ylim = range(-1200,500, finite = TRUE),col='green',asp=asp_ratio,labels =lvls,
        vfont = c("sans serif", "bold"), labcex = 2.5,cex.lab = 2.5,cex.axis = 3.0,lwd = 2, lty = 1, xlab = "X(m)", ylab = "Z(m)")


#legend(x='topright',legend="",
#       col=c("green","red"), pch=c(NA,NA),lty = c(2),  cex=0.8,lwd = 2)
#mtext(paste("Days= 100"), side = 1, line = -16, adj = 0.1, cex = 1.0)
#mtext(paste("Pe=1"), side = 1, line = -14, adj = 0.1, cex = 0.8)   
mtext(paste("(b) XZ view through row",irow), side = 1, line = -45, adj = 0.1, cex = 2.0)
mtext(paste("Vertical Exaggeration :",asp_ratio), side = 1, line = -2.0, adj = 0.9, cex = 1.5)

#YZ
#lvls <-c( 10,20,30)

asp_ratio=2
y<-YZ_ConcLag$YCOORD
z<-YZ_ConcLag$ZCOORD
my.matrix  <- interp(y,z,YZConc)
contour(my.matrix, levels=lvls,xlim = range(4427000:4431000, finite = TRUE),ylim = range(-1200,500, finite = TRUE),col='green',asp=asp_ratio,
        vfont = c("sans serif", "bold"), labcex = 2.5,cex.lab = 2.5,lwd = 2,cex.axis = 3.0, lty = 1, xlab = "Y(m)", ylab = "Z(m)")

# Add single star at specific coordinate
points( 4429475,-50, pch=17, col="red", cex=3.0)

# Add the symbol for the legend
points(4428300.0,-1340, pch=17, col="red", cex=3.0)

#legend(x='topright',legend="",
#       col=c("green","red"), pch=c(NA,NA),lty = c(2),  cex=0.8,lwd = 2)
#mtext(paste("Days= 100"), side = 1, line = -16, adj = 0.1, cex = 1.0)
#mtext(paste("Pe=1"), side = 1, line = -14, adj = 0.1, cex = 0.8)   
mtext(paste("(c) YZ view through column",icol), side = 1, line = -45, adj = 0.1, cex = 2.0)
mtext(paste("Vertical Exaggeration:",asp_ratio ), side = 1, line = -2.0, adj = 0.9, cex = 1.5)
mtext(paste("Release location"), side = 1, line = -2.0, adj = 0.05, cex = 1.5)
dev.off() 

#--End of Contour Plots---------------------------------------------------------------------------------------------------------------------------------------------------!






