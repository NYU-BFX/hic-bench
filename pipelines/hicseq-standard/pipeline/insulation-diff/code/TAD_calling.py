import rpy2.robjects as robjects
import warnings
from rpy2.rinterface import RRuntimeWarning
import numpy as np
import os
import cv2
import math


def TAD_calling(file, up, down):
    warnings.filterwarnings("ignore", category=RRuntimeWarning)
    r = robjects.r
    r('''
        options(warn=-1)
        library("changepoint")
    
        edge_row <- function(l){
            m.binseg <- cpt.meanvar(l,test.stat = "Normal", method = "BinSeg")
            change = cpts(m.binseg)
            if(length(change)>2){
                change = change
            }
            return(change)
        }
    
        edge_detection <- function(y,fold){
            options(warn=-1)
            require(bcp)
            if(fold <= 1){
            p_lim = 0.6
            }
            else{p_lim = 0.4}
            freq = data.frame(table(unlist(y)))
            freq[,1] <- as.numeric(as.character(freq[,1]))
            freq[,2] <- as.numeric(as.character(freq[,2]))
            freq[,3] <- freq[,2]/max(freq[,2])
    
            ff = floor(length(freq[,3])/5000)
            point = freq[1,]
            i = 1
            for(i in 1:(ff+1)){
            f0 = freq[(1+(i-1)*5000):min((i*5000),length(freq[,3])),]
            bcp.1 <- bcp(f0[,3],mcmc = 5000)
            pre.p = data.frame(bcp.1$posterior.prob, bcp.1$posterior.mean[,1])                       
            colnames(pre.p) = c('Probability', 'X1')
            point = as.matrix(f0[pre.p[,1]>0.2,])
            point = point[!is.na(point[,1]),]
            }
            if(!is.null(nrow(point))){
            return(point[-1,])
            }
            else{
            return(point)
            }
        }
    
        boundary_probability <- function(points,cutoff){
          options(warn=-1)
          require('cluster')
          points = as.data.frame(points)
          d=dist(points[,1],method = "euclidean")
          if(nrow(points)>0){
            cluster = hclust(d,method = "complete")
            clusterCut <- cutree(cluster, h = cutoff)
            points[,4] = clusterCut
            pk1 = c()
            j=1
            for(i in unique(points[,4])){
              m = max(points[points[,4]==i,2])
              pk1[j]=points[points[,2]==m&points[,4]==i,1]
              j=j+1
            }
            return(pk1)
          }
          else{
            return(points)
          }
        }
    
        boundary_detection = function(file_path,cover_lim_up,cover_lim_down,fold) {
          if(fold<=1){
            range = 500
            range2 = 150
            range3 = 20
          }
          else{
            range = 500/fold
            range2 = 150/fold
            range3 = 20/fold
          }
    
          bps = data.frame(0,0,0,0,0,0)
          contact_map = file(file_path, "r")
          line = readLines(contact_map, n = 1)
          l = as.numeric(as.character(strsplit(line,'\t')[[1]]))
          end = length(l) - 50
          i = 1
          while ( i<end ) {
            line = readLines(contact_map, n = 1)
            l = as.numeric(as.character(strsplit(line,'\t')[[1]]))
            l1=l[i:min((i+range),length(l))]
            l1[l1>cover_lim_up]=cover_lim_up
            l1[l1<=cover_lim_down]=0
            l2=l[i:max((i-range),1)]
            l2[l2>cover_lim_up]=cover_lim_up
            l2[l2<=cover_lim_down]=0
            if ( length(line) == 0 ) {
              break
            }
            else{
              if(i>50){
                edge_right = edge_row(l1)
                edge_left = edge_row(l2)
                if(length(edge_right)>0&length(edge_left)>0){
                  bps[i,1] = i
                  bps[i,2] = i + edge_right[2]
                  bps[i,3] = i - edge_left[2]
                  bps[i,4] = i + edge_right[3]
                  bps[i,5] = i - edge_left[3]
                }
              }
            }
            i = i+1
          }
          bps[,1] = 1:length(bps[,1])
          boundary_col = bps[(bps[,2]-bps[,1])<range2&(bps[,2]-bps[,1])>range3,2]
          boundary_col = boundary_col[!is.na(boundary_col)]
          boundary_row = bps[(bps[,1]-bps[,3])<range2&(bps[,1]-bps[,3])>range3,3]
          boundary_row = boundary_row[!is.na(boundary_row)]
          close(contact_map)
          out_list = list(bps,boundary_col,boundary_row)
          return(out_list)
        }
    
        TAD_detection<-function(boundary_left,boundary_right,band1,distance_cutoff,fold){
          options(warn=-1)
          require(MASS) 
          require(spatstat)
          prob_left = edge_detection(boundary_left,fold)
          if(FALSE){
            return(data.frame(0,0,0))
          }
          else{
            if(is.null(nrow(prob_left)) | nrow(prob_left)==0){
              return(data.frame(0,0,0))
            }
            else{
              prob_left2 = prob_left[order(prob_left[,1]),]
              edg_col = boundary_probability(prob_left2,distance_cutoff)
              prob_right = edge_detection(boundary_right,fold)
              if(!is.null(nrow(prob_right))){
                if(nrow(prob_right)>0){
                  prob_right2 = prob_right[order(prob_right[,1]),]
                  edg_row = boundary_probability(prob_right2,distance_cutoff)
                  pnts = data.frame(0,0)
                  t=1
                  for(i in 1:length(edg_col)){
                    for(j in 1:length(edg_row)){
                      if(edg_col[i]-edg_row[j]>10&edg_col[i]-edg_row[j]<1000){
                        pnts[t,1]=edg_col[i]
                        pnts[t,2]=edg_row[j]
                        t=t+1
                      }
                    }
                  }
    
                  coordinate = band1
                  coordinate = na.omit(coordinate)
                  pre_boundary=data.frame()
    
                  for(t in 1:ceiling(length(coordinate[,1])/3000)){
                    n1=1+3000*(t-1)
                    n2=3000*t
                    x = na.omit(coordinate[n1:min(n2,length(coordinate[,1])),])
                    xy.kde <- kde2d(x[,1],x[,2], n=(n2-n1)/3,h=20)
                    xy.im <- im(t(xy.kde$z), xcol=xy.kde$x, yrow=xy.kde$y)
                    pnt_tmp = pnts[pnts[,2]<xy.im$xrange[2]&pnts[,2]>xy.im$xrange[1],]
                    pnt_tmp = pnt_tmp[pnt_tmp[,1]<xy.im$yrange[2]&pnt_tmp[,1]>xy.im$yrange[1],]
                    c.0 <- sum(xy.kde$z)
                    pre_boundary_temp = data.frame()
                    if(length(pnt_tmp[,1])>0){
                      for(i in 1:length(pnt_tmp[,1])){
                        z <- interp.im(xy.im, pnt_tmp[i,2], pnt_tmp[i,1])                            
                        p <- sapply(z, function(a) sum(xy.kde$z[xy.kde$z < a])) / c.0
                        pre_boundary_temp[i,1] = p
                        pre_boundary_temp[i,2] = pnt_tmp[i,2]
                        pre_boundary_temp[i,3] = pnt_tmp[i,1]
                      }
                      pre_boundary = rbind(pre_boundary,pre_boundary_temp)
                    }
                  }
                  if(fold<=1){
                    boundary = pre_boundary[pre_boundary[,1]>0.8,]
                  }
                  else(boundary = pre_boundary[pre_boundary[,1]>quantile(pre_boundary[,1],0.8),])
                  return(boundary)
                }
                else{
                  return(data.frame(0,0,0))
                }
              }
              else{
                return(data.frame(0,0,0))
              }
            }
          }
        }
    
        TAD_caller <- function(file, color_lim_up, color_lim_down, distance_cutoff = 25, fold = 1){
          options(warn=-1)
          require(MASS)
          require(spatstat) 
          boundary_list = boundary_detection(file, color_lim_up,color_lim_down, fold)
          boundary_left = boundary_list[[2]]
          boundary_right = boundary_list[[3]]
          if(length(boundary_right)==0|length(boundary_left)==0){
            return(t(matrix(c(0,0,0))))
          }
          else{
            pair_boundary1 = boundary_list[[1]][,c(3,2)]
            pair_boundary2 = boundary_list[[1]][,c(5,4)]
            boundary = t(matrix(c(1,1,1)))
            for(i in 1:max(1,floor(length(boundary_left)/2000))){
              b1 = boundary_left[(1+(i-1)*2000):min((i*2000),length(boundary_left))]
              b2 = boundary_right[(1+(i-1)*2000):min((i*2000),length(boundary_left))]
              b1 = na.omit(b1)
              b2 = na.omit(b2)
              b3 = pair_boundary1[min(b1,b2):max(b1,b2),]
              b4 = pair_boundary2[min(b1,b2):max(b1,b2),]
              b5 = rbind(as.matrix(b3),as.matrix(b3),as.matrix(b4))
              b5 = b5[order(b5[,1]),]
              boundary_temp = TAD_detection(b1,b2,b5,distance_cutoff,fold)
              boundary = rbind(boundary,as.matrix(boundary_temp))
            }
            return(boundary[-1,])
          }
        }
        ''')
    TAD_calling = r['TAD_caller']
    try:
        contact_map = np.loadtxt(file)
        if np.sum(contact_map) > 0:
            TAD = np.array(TAD_calling(file, up, down))
        else:
            TAD = np.array([0, 0, 0])
    except:
        TAD = np.array([0, 0, 0])
    return(TAD)


def file_split(file, outfile, length=1000, print_subcontact=0, aliases='None'):
    '''
    divide whole chromosome into different sub chromosome files and calculate cutoffs
    :param file: input chromosome file
    :param outfile: output path for sub chromosome files
    :param length: length of sub chromosome
    :param print_subcontact: if print_subcontact==1, not save sub chromosome files, if print_subcontact!=1, save sub
    chromosome files
    :return: up limit and down limit for this chromosome
    '''
    base = os.path.basename(file)
    if aliases != "None":
        filename = aliases
    else:
        filename = os.path.splitext(base)[0]

    print(file)
    contact_map = np.loadtxt(file)
    r, c = contact_map.shape
    num = int(r/length)
    mean_list = []
    files = []

    flag = 0
    if (num*length+length/2) > c:
        flag = 1

    if flag == 1:
        for i in range(num):
            sub_contact = contact_map[(i * length):((i + 1) * length), (i * length):((i + 1) * length)]
            out = os.path.join(outfile, filename + '.' + str(i * length) + '.' + str((i + 1) * length) + '.subchr')
            files.append(out)
            mean_list.append(np.mean(sub_contact))
            if print_subcontact != 1:
                np.savetxt(out, sub_contact, delimiter='\t')
            if i != (num-1):
                sub_contact2 = contact_map[int(i * length+length/2):int((i + 1) * length+length/2),
                               int(i * length+length/2):int((i + 1) * length+length/2)]
                out2 = os.path.join(outfile,
                                   filename + '.' + str(int(i * length+length/2)) + '.' + str(int((i + 1) * length + length/2)) + '.subchr')
                files.append(out2)
                mean_list.append(np.mean(sub_contact2))
                if print_subcontact != 1:
                    np.savetxt(out2, sub_contact2, delimiter='\t')
        sub_contact = contact_map[(num * length):c, (num * length):c]
        out = os.path.join(outfile, filename + '.' + str(num * length) + '.' + str(c) + '.subchr')
        files.append(out)
        mean_list.append(np.mean(sub_contact))
        if print_subcontact != 1:
            np.savetxt(out, sub_contact, delimiter='\t')
        sub_contact2 = contact_map[int((num-1) * length + length / 2):int(((num-1) + 1) * length + length / 2),
                       int((num-1) * length + length / 2):int(((num-1) + 1) * length + length / 2)]
        out2 = os.path.join(outfile, filename + '.' + str(int((num-1) * length + length / 2)) + '.' + str(int(((num-1) + 1) * length + length / 2)) + '.subchr')
        files.append(out2)
        mean_list.append(np.mean(sub_contact2))
        if print_subcontact != 1:
            np.savetxt(out2, sub_contact2, delimiter='\t')
    else:
        for i in range(num):
            sub_contact = contact_map[(i * length):((i + 1) * length), (i * length):((i + 1) * length)]
            out = os.path.join(outfile, filename + '.' + str(i * length) + '.' + str((i + 1) * length) + '.subchr')
            files.append(out)
            mean_list.append(np.mean(sub_contact))
            if print_subcontact != 1:
                np.savetxt(out, sub_contact, delimiter='\t')
            sub_contact2 = contact_map[int(i * length+length/2):int((i + 1) * length+length/2),
                           int(i * length+length/2):int((i + 1) * length+length/2)]
            out2 = os.path.join(outfile,
                               filename + '.' + str(int(i * length+length/2)) + '.' + str(int((i + 1) * length + length/2)) + '.subchr')
            files.append(out2)
            mean_list.append(np.mean(sub_contact2))
            if print_subcontact != 1:
                np.savetxt(out2, sub_contact2, delimiter='\t')
        sub_contact = contact_map[(num * length):c, (num * length):c]
        out = os.path.join(outfile, filename + '.' + str(num * length) + '.' + str(c) + '.subchr')
        files.append(out)
        mean_list.append(np.mean(sub_contact))
        if print_subcontact != 1:
            np.savetxt(out, sub_contact, delimiter='\t')
        sub_contact2 = contact_map[int(num * length + length / 2):int((num + 1) * length + length / 2),
                       int(num * length + length / 2):int((num + 1) * length + length / 2)]
        out2 = os.path.join(outfile,
                            filename + '.' + str(int(num * length + length / 2)) + '.' + str(int(
                                (num + 1) * length + length / 2)) + '.subchr')
        files.append(out2)
        mean_list.append(np.mean(sub_contact2))
        if print_subcontact != 1:
            np.savetxt(out2, sub_contact2, delimiter='\t')
    mean_list = np.array(mean_list)
    mean_list = mean_list[~np.isnan(mean_list)]
    cutoff_up = np.median(np.array(mean_list) / 0.22)
    cutoff_down = np.median(np.array(mean_list) / 2)
    files.sort()
    return cutoff_up, cutoff_down, files


def TAD_plot(f1, output, TAD, down, up):
    '''
    :param f1: contact map
    :param down: down limit for color key
    :param up: up limit for color key
    :return:
    '''
    img = np.loadtxt(f1)
    img[img > up] = up
    img[img < down] = 0
    img = img * (255 / up)
    img = img.astype(np.uint8)
    color_img = cv2.cvtColor(img, cv2.COLOR_GRAY2BGR)
    color_img[:, :, 0] = 255 - color_img[:, :, 0]
    color_img[:, :, 1] = 255 - color_img[:, :, 1]
    color_img[:, :, 2] = 255
    for x in TAD:
        if float(x[0]) > 0.8:
            color_img = cv2.rectangle(color_img, (math.ceil(float(x[1])), math.ceil(float(x[1]))),
                                      (math.ceil(float(x[2])), math.ceil(float(x[2]))), (0, 0, 0), 2,
                                      lineType=4)
    cv2.imwrite(output, color_img)


def subTAD_calling(path, up, down, plt=1, prnt=1, aliase='None'):
    dirs = os.listdir(path)
    tad = np.empty(shape=[0, 3])
    for file in dirs:
        if file[-7:] == '.subchr':

            print(file)

            start = file.split('.')[-3]
            tad_temp = TAD_calling(os.path.join(path, file), up, down)
            tad_location = tad_temp

            if plt == 1:
                try:
                    TAD_plot(os.path.join(path, file), os.path.join(path, file + '.tiff'),
                                         tad_temp, down, up)
                except:
                    pass
            try:
                tad_location[:, 1:3] = tad_location[:, 1:3] + int(start)
                if not np.array_equal(tad_location, np.array([0, int(start), int(start)])):
                    tad = np.vstack((tad, tad_location))
                if prnt == 1:
                    np.savetxt(os.path.join(path, file + '.TAD.txt'), tad_location, delimiter="\t")
            except:
                tad_location = np.array([[0, 0, 0]])
                if prnt == 1:
                    np.savetxt(os.path.join(path, file + '.TAD.txt'), tad_location, delimiter="\t")
    sorted_idx = np.lexsort(tad[:, 1:3].T)
    sorted_tad = tad[sorted_idx, :]
    row_mask = np.append([True], np.any(np.diff(sorted_tad[:, 1:3], axis=0), 1))
    unique_tad = sorted_tad[row_mask]
    if aliase == 'None':
        np.savetxt(os.path.join(path, 'TAD.merge.txt'), unique_tad, delimiter="\t")
    else:
        np.savetxt(os.path.join(path, aliase + '.TAD.merge.txt'), unique_tad, delimiter="\t")
