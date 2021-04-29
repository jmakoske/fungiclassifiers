#!/usr/bin/env python3
# FILE: trainCNN.py
# AUTHOR: Duong Vu
# CREATE DATE: 07 June 2019
import sys
#if sys.version_info[0] >= 3:
#        unicode = str
import os, argparse
import tensorflow as tf
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, InputLayer
from tensorflow.keras.layers import Convolution1D, MaxPooling1D
from tensorflow.keras.layers import Dropout, Activation, Flatten
from tensorflow.keras import utils
from tensorflow.keras import backend as K
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score as acc
from sklearn.metrics import matthews_corrcoef as mcc
import numpy as np
import json
import time

parser=argparse.ArgumentParser(prog='trainCNN.py', 
                                                           usage="%(prog)s [options] -i fastafile -c classificationfile -p classificationposition",
                                                           description='''Script that trains a CNN model to classify sequences''',
                                                           epilog="""Written by Duong Vu duong.t.vu@gmail.com""",
   )

parser.add_argument('-i','--input', required=True, help='the fasta file')
parser.add_argument('-o','--out', help='The folder name containing the model and associated files.') #optional
parser.add_argument('-c','--classification', required=True, help='the classification file in tab. format.')
parser.add_argument('-p','--classificationpos', required=True, type=int, default=0, help='the classification position to load the classification.')
parser.add_argument('-k','--kmer', type=int, default=6, help='the k-mer for the representation of the sequences.')

args=parser.parse_args()
fastafilename= args.input
classificationfilename=args.classification
classificationlevel=args.classificationpos
k = args.kmer
modelname=args.out

#fastafilename=sys.argv[1]
#classificationfilename=sys.argv[2] #taxonomy file 
#classificationlevel=int(sys.argv[3]) #the level of classification to get taxa from the taxonomy file
#k = 6
#if len(sys.argv) >4:
#       k= int(sys.argv[4])

def GetBase(filename):
        return filename[:-(len(filename)-filename.rindex("."))]
                
def load_data(fastafilename,classificationfilename,classificationlevel):
        #load classification
        #allseqids=[]
        print('  Reading', classificationfilename)
        records= open(classificationfilename, encoding = "ISO-8859-1")
        classification=[]
        classificationdict={}
        level=""
        for record in records:
                texts=record.split("\t")
                if record.startswith("#"):
                        if classificationlevel < len(texts):
                                level=texts[classificationlevel].rstrip()
                        continue                
                seqid=texts[0].replace(">","").rstrip()
                classname=""
                if classificationlevel < len(texts):
                        classname=texts[classificationlevel].rstrip()
                if classname !="":
                        #allseqids.append(seqid)
                        classification.append(classname)
                        classificationdict[seqid]=classname
        records.close()
        print('  Done')#, len(allseqids))
        classificationset=set(classification)
        classes=list(classificationset)
        #load fastafile, save a new fasta file containing only sequences having a classification
        fastafile=open(fastafilename)
        newfastafilename=GetBase(fastafilename) + "." + str(classificationlevel) + ".fasta"
        newfastafile=open(newfastafilename,"w")
        writetofile=False
        seq=""
        sequences=[]
        seqids=[]
        taxa=[]
        tic = time.time()
        for i_l, line in enumerate(fastafile):
                if i_l%5000==0:
                        tac = time.time()
                        print(i_l, tac-tic)
                        tic = tac
                if line.startswith(">"):
                        if writetofile==True:
                                sequences.append(seq)
                        writetofile=False
                        seq=""
                        seqid=line.split("|")[0].replace(">","").rstrip()
                        #print(seqid)
                        if seqid in classificationdict: 
                                taxonname=classificationdict[seqid]
                                if taxonname !="":
                                        newfastafile.write(">" + seqid + "\n")
                                        seqids.append(seqid)
                                        taxa.append(taxonname)
                                        writetofile=True
                else:
                        if writetofile==True:
                                newfastafile.write(line)
                                seq=seq+line.rstrip()
        if writetofile==True:
                sequences.append(seq)
        fastafile.close()
        newfastafile.close()
        return newfastafilename, seqids, sequences, taxa, classes,level

def load_matrix(matrixfilename,seqids,sequences,taxa,classes):
        #load matrix
        vectors= list(open(matrixfilename, "r"))
        vectors=vectors[1:]
        print(len(vectors), 'vectors')
        X=[]
        Y=[]
        seqIDList=[]
        seqList=[]
        seqidsdict = {}
        for i_s, s in enumerate(seqids):
                seqidsdict[s]=i_s
        for i in range(0,len(classes)):
                seqIDList.append([])
                seqList.append([])
        tic = time.time()
        for i_v, vector in enumerate(vectors):
                if i_v%5000==0:
                        tac = time.time()
                        print(i_v, tac-tic)
                        tic = tac
                #print(vector)
                elements=vector.split(",")
                seqid=elements[0]
                #taxonname= taxa[seqids.index(seqid)]
                taxonname = taxa[seqidsdict[seqid]]
                #seq=sequences[seqids.index(seqid)]
                seq = sequences[seqidsdict[seqid]]
                index=classes.index(taxonname)
                #print(seqid, taxonname, seq, index)
                Y.append(index)
                X.append(elements[1:])
                #if seqid not in seqIDList[index]:
                        #seqIDList[index].append(seqid)
                        #seqList[index].append(seq)
        X=np.array(X,dtype=float)
        Y=np.array(Y,dtype=int)
        data_max=0
        #data_max= np.amax(X)
        #X = X/data_max
        return X,Y,len(classes),len(X[0]),data_max,seqIDList,seqList

def multi_mcc(y_true, y_pred):
        confusion_m = tf.matmul(K.transpose(y_true), y_pred)

        N = K.sum(confusion_m)

        up = N*tf.linalg.trace(confusion_m) - K.sum(tf.matmul(confusion_m, confusion_m))
        down_left = K.sqrt(N**2 - K.sum(tf.matmul(confusion_m, K.transpose(confusion_m))))
        down_right = K.sqrt(N**2 - K.sum(tf.matmul(K.transpose(confusion_m), confusion_m)))
        
        mcc = up / (down_left * down_right + K.epsilon())
        mcc = tf.where(tf.math.is_nan(mcc), tf.zeros_like(mcc), mcc)

        return K.mean(mcc)
                
def create_model(nb_classes,input_length):
        model = Sequential()
        print('input_length=', input_length)
        model.add(InputLayer(input_shape=(input_length, 1)))
        model.add(Convolution1D(5,5, padding='valid')) #input_dim
        model.add(Activation('relu'))
        model.add(MaxPooling1D(pool_size=2,padding='valid'))
        model.add(Convolution1D(10, 5,padding='valid'))
        model.add(Activation('relu'))
        model.add(MaxPooling1D(pool_size=2,padding='valid'))
        model.add(Flatten())
        ##
        ##MLP
        model.add(Dense(500))
        model.add(Activation('relu'))
        model.add(Dropout(0.5))
        model.add(Dense(nb_classes))
        model.add(Activation('softmax'))
        model.compile(optimizer='adam', loss='categorical_crossentropy',
                      metrics=['accuracy', multi_mcc])
        print(model.summary())
        return model

def SaveConfig(configfilename,classifiername,fastafilename,jsonfilename,classificationfilename,classificationpos,kmer,data_max):
        if not classifiername.startswith("/"):
                classifiername=os.getcwd() + "/" + classifiername
        if not fastafilename.startswith("/"):
                fastafilename=os.getcwd() + "/" + fastafilename
        if not jsonfilename.startswith("/"):
                jsonfilename=os.getcwd() + "/" + jsonfilename
        if not classificationfilename.startswith("/"):
                classificationfilename=os.getcwd() + "/" + classificationfilename
        model="cnn"
        #save the config:classifierfilename, model, classificationfilename,classificationpos,k-mer
        configfile=open(configfilename,"w")
        configfile.write("Classifier name: " + classifiername + "\n")
        configfile.write("Model: " + model + "\n")
        configfile.write("Data max: " + str(data_max) + "\n")
        configfile.write("K-mer number: " + str(kmer) + "\n")
        configfile.write("Fasta filename: " + fastafilename + "\n")
        configfile.write("Classification filename: " + classificationfilename + "\n")
        configfile.write("Column number to be classified: " + str(classificationpos) + "\n")
        configfile.write("Classes filename: " + jsonfilename + "\n")
        configfile.close()      

def SaveClasses(jsonfilename,classnames,seqIDList,seqList):
        #create json dict
        taxadict={}
        i=0
        count=0
        for classname in classnames:    
                #classname=unicode(classname,errors='ignore')
                currentclass=[]
                j=0        
                for seqID in seqIDList[i]:
                        seq=seqList[i][j]
                        currentclass.append({'id': seqID, 'seq': seq})
                        j=j+1
                        count=count+1            
                taxadict[classname]=currentclass
                i=i+1
        #write to file
        with open(jsonfilename,"w") as json_file:
                json.dump(taxadict,json_file) #,encoding='latin1')

if __name__ == "__main__":
        path=sys.argv[0]
        path=path[:-(len(path)-path.rindex("/")-1)]
        #load data
        print('Loading FASTA data from', fastafilename, 'and')
        print('classification data from', classificationfilename, classificationlevel)
        newfastafilename,seqids,sequences,taxa,classes,level = load_data(fastafilename,classificationfilename,classificationlevel)
        #represent sequences as matrix of k-mer frequencies
        filename=GetBase(fastafilename)
        matrixfilename=filename + "." + str(k) + ".matrix"
        command=path + "fasta2matrix.py " +  str(k) + " " + newfastafilename + " " + matrixfilename
        print('Running command:', command)
        os.system(command)
        print('Loading matrix data from', matrixfilename)
        X,Y,nb_classes,input_length,data_max,seqIDList,seqList = load_matrix(matrixfilename,seqids,sequences,taxa,classes)
        #sys.exit()
        #train data
        trainids=np.array(range(0,len(X)),dtype=int)
        traindata=X[trainids]
        trainlabels=Y[trainids]
        #training
        print('Creating data with nb_classes={}, input_length={}'.format(nb_classes, input_length))
        model = create_model(nb_classes,input_length)
        traindata = traindata.reshape(traindata.shape + (1,))
        trainlabels_bin=utils.to_categorical(trainlabels, nb_classes)
        print('Training data: X:', traindata.shape, 'Y:', trainlabels_bin.shape, 'actual classes:', len(np.unique(trainlabels)))
        #model.fit(traindata, trainlabels_bin, validation_split=0.2,
        #          epochs=10, batch_size=20, verbose = 2)
        x_train, x_valid, y_train, y_valid = train_test_split(traindata, trainlabels_bin,
                                                              test_size=0.2, shuffle=True)
        model.fit(x_train, y_train, validation_data=(x_valid, y_valid),
                  epochs=10, batch_size=20, verbose=2)
        pred_train = np.argmax(model.predict(x_train), axis=1)
        pred_valid = np.argmax(model.predict(x_valid), axis=1)
        y_train = np.argmax(y_train, axis=1)
        y_valid = np.argmax(y_valid, axis=1)
        #print(y_train)
        #print(pred_train)
        print('Train: accuracy={}, mcc={}'.format(acc(y_train, pred_train), mcc(y_train, pred_train)))
        print('Valid: accuracy={}, mcc={}'.format(acc(y_valid, pred_valid), mcc(y_valid, pred_valid)))
        #save model
#       modelname=filename.replace(".","_") + "_cnn_classifier"
        if modelname==None or modelname=="":
                modelname=filename.replace(".","_") + "_cnn_classifier" 
                if level !="":
                        modelname=filename + "_" + level + "_cnn_classifier"
        basename=modelname
        if "/" in modelname:
                basename=modelname[modelname.rindex("/")+1:]
        if os.path.isdir(modelname) == False:
                os.system("mkdir " + modelname)
        #save model
        classifiername=modelname + "/" + basename + ".classifier"
        model.save(classifiername)
        #save seqids for each classification
        jsonfilename=modelname + "/" + basename + ".classes"
        SaveClasses(jsonfilename,classes,seqIDList,seqList)
        #save config    
        configfilename=modelname + "/" + basename + ".config"
        SaveConfig(configfilename,classifiername,fastafilename,jsonfilename,classificationfilename,classificationlevel,k,data_max)
        print("The classifier is saved in the folder " + modelname + ".")
        
