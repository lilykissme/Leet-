from mesh_to_mesh_trie import *
import pandas as pd
import numpy as np
from collections import Counter
import logging
import time
import multiprocessing


def getMeshNumber(Meshtree, disName):
	disName = disName
	meshDict = Meshtree.string_to_mesh_dictionary(disName)
	semSim = []
	point = []
	if len(meshDict) == 0:
		# print("disease not in Mesh")
		pass
	else:
		disName = list(meshDict.keys())
		meshNumber = meshDict[disName[0]]

		for i in range(len(meshNumber)):
			temp = meshNumber[i].split('.')
			for j in range(len(temp)):
				semSim.append(pow(0.5, j))

				point.append('.'.join(temp[:len(temp) - j]))
	return semSim, point


def letBegin():
	Mesh = MeSHTrie()
	allDisList = np.loadtxt(r'..\..\hmddV3\alldisList.txt', delimiter='\n', dtype=np.str_)

	newDisList = []

	for dis in allDisList:
		meshDict = Mesh.string_to_mesh_dictionary(dis)
		if len(meshDict) == 0:
			pass
		else:
			disName = list(meshDict.keys())[0]
			if disName not in newDisList:
				newDisList.append(disName)

	# np.savetxt(r'.\source\hmddV3\newdis_Name.txt',newDisList,fmt='%s', delimiter='\n',encoding='utf-8')

	disList = np.array(newDisList)

	print(len(Mesh.mesh_dict.values()))
	'''
	for i in Mesh.mesh_dict.keys():
		print(Mesh.mesh_dict[i],i)
		break

	print(Mesh.string_to_mesh_dictionary(disList[1]))
	print(Mesh.string_to_mesh_dictionary("acute heart failure"))
	'''

	# TODO
	# newDisList--- 计算相似性矩阵。
	# in Mesh tree  and newDisList...

	allData = pd.read_excel(r'..\..\alldata.xlsx')

	allmirList = allData.sort_values(by='mir').loc[:, 'mir']
	# alldisList=allData.sort_values(by='disease').loc[:,'disease']
	mdData = allData.loc[:, ['mir', 'disease']]

	mirList = []
	# disList=[]
	count = 0
	for item in allmirList:
		item = item.lower()
		if item not in mirList:
			mirList.append(item)
		count += 1

	# dis-gen data
	disGen = pd.read_csv(r'..\..\disGen\all_gene_disease_associations.tsv', sep='\t', header=0)

	dis_G = disGen.loc[:, ['diseaseName']]
	gen_G = disGen.loc[:, ['geneId']]
	assocation = disGen.loc[:, ['geneId', 'diseaseName']]

	disList_G = []
	genList_G = []
	notIn = []

	count_dis = 0
	count_gen = 0
	num = assocation.index.tolist()

	for i in num:
		itemD = dis_G.iloc[i, 0]

		meshDict = Mesh.string_to_mesh_dictionary(itemD)
		if len(meshDict) == 0:
			pass
		else:
			itemD = list(meshDict.keys())[0]

			if (itemD not in disList_G):
				disList_G.append(itemD)
				count_dis += 1

		itemG = gen_G.iloc[i, 0]

		if itemG not in genList_G:
			genList_G.append(itemG)
			count_gen += 1

	right_D = 0
	for i in range(len(disList)):

		temp1 = disList[i].lower()

		for j in range(len(disList_G)):
			temp2 = disList_G[j].lower()
			disList_G[j] = disList_G[j].lower()
			if temp1 == temp2:
				# if temp in disList_G:
				notIn.append(disList_G[j])

				right_D += 1

	notInN = np.array(notIn)
	np.savetxt(r'..\..\hmddV3\NewdisList.txt', notInN, fmt='%s', delimiter='\n')

	mdMatrix = np.zeros((len(mirList), len(notInN)))

	pos = []
	record = []
	right = []
	for i in range(count):
		mir = mdData.iloc[i, 0].lower()
		dis = mdData.iloc[i, 1].lower()

		meshDict = Mesh.string_to_mesh_dictionary(dis)
		if len(meshDict) == 0:
			pass
		else:
			dis = list(meshDict.keys())[0]

			if dis in notIn:
				col = mirList.index(mir)
				row = notIn.index(dis)
				pos.append(col * len(mirList) + row)

				mdMatrix[col, row] = 1

	mirListN = np.array(mirList)

	print("semsim!")
	disSemsim = np.zeros((len(notInN), len(notInN)))

	for disA in notInN:
		for disB in notInN:
			col = int(np.where(notInN == disA)[0])
			row = int(np.where(notInN == disB)[0])
			if disA == disB:
				disSemsim[col, row] = 1
			else:
				semSimA, pointA = getMeshNumber(Mesh, disA)
				semSimB, pointB = getMeshNumber(Mesh, disB)

				samePoint = []
				if len(semSimA) != 0 and len(semSimB) != 0:
					for ancestorA in pointA:
						for ancestorB in pointB:
							if ancestorB == ancestorA:
								indexA = pointA.index(ancestorA)
								indexB = pointB.index(ancestorB)
								samePoint.append(semSimA[indexA])
								samePoint.append(semSimB[indexB])

					disSemsim[col, row] = sum(samePoint) / (sum(semSimA) + sum(semSimB))

	np.savetxt(r'..\..\hmddV3\disSemsim.txt', disSemsim, fmt='%f', delimiter='\t', encoding='utf-8')

	np.savetxt(r'..\..\hmddV3\NewallData.txt', mdMatrix, fmt='%d', delimiter='\t')

	dis_Gen = {}
	# dis_genIndex
	dis_GenIndex = {}

	genList_G.sort()
	genList_GNd = np.array(genList_G)

	for i in notInN:
		dis_Gen[i] = []
		dis_GenIndex[i] = []

	genDisNum = []
	for i in num:
		dis = assocation.iloc[i, 1].lower()
		gen = assocation.iloc[i, 0]

		meshDict = Mesh.string_to_mesh_dictionary(dis)
		if len(meshDict) == 0:
			pass
		else:
			dis = list(meshDict.keys())[0]
			if dis in dis_Gen.keys() and gen not in dis_Gen[dis]:
				dis_Gen[dis].append(gen)

	for key in dis_Gen.keys():
		genDisNum.append(str(len(dis_Gen[key])) + "-" + key)


	gLLS = pd.read_csv(r'..\..\HumanNet.v1.join.txt', sep='\t', header=None)

	# print(gLLS.head())
	gen1 = gLLS.iloc[:, 0]
	gen2 = gLLS.iloc[:, 1]
	gLLSNum = int(gLLS.shape[0])
	gen_Gen = {}

	gen1Num = 0
	gen2Num = 0
	gen1Name = []
	gen2Name = []

	for i in range(gLLSNum):
		temp1 = int(gen1.iloc[i])
		temp2 = int(gen2.iloc[i])
		if temp1 not in gen1Name:
			gen1Name.append(temp1)
			gen1Num += 1
		if temp2 not in gen2Name:
			gen2Name.append(temp2)
			gen2Num += 1

	gen1Name.sort()
	gen2Name.sort()
	gLLsMat = np.zeros((gen1Num, gen2Num))

	conbamGen = []
	for i in gen1Name:
		if i in genList_GNd:
			if i not in conbamGen:
				conbamGen.append(i)

	for i in gen2Name:
		if i in genList_GNd:
			if i not in conbamGen:
				conbamGen.append(i)
	conbamGen.sort()
	print("gen complete!")
	lenGen = len(conbamGen)
	genList_GMat = np.zeros((lenGen, lenGen))

	conbamGenNd = np.array(conbamGen)

	for i in range(len(gen1Name)):
		# gen_Gen[gen1Name[i]]={}
		gen_Gen[gen1Name[i]] = []

	for i in range(gLLSNum):
		temp1 = int(gen1.iloc[i])
		temp2 = int(gen2.iloc[i])
		tempLLS = float(gLLS.iloc[i, 23])

		gLLsMat[gen1Name.index(temp1), gen2Name.index(temp2)] = tempLLS

		if temp1 in gen_Gen.keys():
			gen_Gen[temp1].append(temp2)

		if temp1 in conbamGenNd and temp2 in conbamGenNd:
			gCol = conbamGen.index(temp1)
			gRow = conbamGen.index(temp2)
			genList_GMat[gCol, gRow] = tempLLS
			genList_GMat[gRow, gCol] = genList_GMat[gCol, gRow]

	maxLLS = 4
	for i in range(lenGen):
		genList_GMat[i, i] = maxLLS + 1

	print(maxLLS)
	for i in dis_Gen.keys():
		tempGen = dis_Gen[i]
		for j in tempGen:
			if j in conbamGenNd:
				tempGenIndex = conbamGen.index(j)
				dis_GenIndex[i].append(tempGenIndex)

	# print(gen1Num,gen2Num,gLLSNum)
	print("ready")
	return notIn, dis_Gen, gen1Name, gen2Name, gLLsMat, gen_Gen, genList_GMat, conbamGen, dis_GenIndex


def getDisFunSimDemo2(dis_Gen, gLLsMat, disList, gen1Name, gen2Name, begin, end, gen_GenDict, genListMat, conbam,
					  dis_GenIndex):
	test1 = time.process_time()
	print("当前子进程的名称为%s" % (multiprocessing.current_process()))
	notIn = disList

	disNum = len(notIn)
	disFunS = np.zeros((disNum, disNum))
	gen1NameNd = np.array(gen1Name)
	gen2NameNd = np.array(gen2Name)

	gen_Dict = gen_GenDict

	conbamNd = conbam
	for i in range(begin, end, 1):
		dis1 = notIn[i]
		dis1Gen = dis_GenIndex[dis1]
		dis1GenNd = np.array(dis1Gen)
		for j in range(i + 1):
			dis2 = notIn[j]
			dis2Gen = dis_GenIndex[dis2]
			dis2GenNd = np.array(dis2Gen)

			if dis1 != dis2:
				totalLLS1 = 0
				totalLLS2 = 0

				m = len(dis1Gen)
				n = len(dis2Gen)
				tempMatrix = genListMat[dis1GenNd, :]
				tempMatrix = tempMatrix[:, dis2GenNd]

				if len(tempMatrix) != 0:
					colMax = tempMatrix[range(tempMatrix.shape[0]), np.argmax(tempMatrix, axis=1)]
					# m=len(np.extract(colMax==0,colMax))
					totalLLS1 = np.sum(colMax)

					rowMax = tempMatrix[np.argmax(tempMatrix, axis=0), range(tempMatrix.shape[1])]
					# n=len(np.extract(rowMax==0,rowMax))
					totalLLS2 = np.sum(rowMax)

				if m != 0 or n != 0:
					disFunS[i, j] = (totalLLS1 + totalLLS2) / (m + n)

			else:
				disFunS[i, j] = 1

			disFunS[j, i] = disFunS[i, j]
			# print(dis1, dis2, disFunS[i, j])

			print("one loop", i, j,  time.process_time() - test1)

		np.savetxt(r'.\source\Matrix\Demo2-disfunS-' + str(i) + r'.txt', disFunS, fmt='%f', delimiter='\t',
				   encoding='utf-8')
		print("complete %d" % i)
	np.savetxt(r'.\source\Matrix\Demo2disfunS-' + str(begin) + '-' + str(end) + r'.txt', disFunS, fmt='%f',
			   delimiter='\t',
			   encoding='utf-8')

	print("%snum is %d to %d" % ('mission', begin, end))


# return disFunS

if __name__ == "__main__":  # win操作系统需要加上,否则会出现异常报错RuntimeError

	print("let start!")
	startTime = time.process_time()
	notIn, dis_Gen, gen1Name, gen2Name, gLLsMat, gen_Gen,genList_GMat,conbamGen,dis_GenIndex = letBegin()
	print(time.process_time() - startTime)
	lenth = len(notIn)
	p1 = multiprocessing.Process(target=getDisFunSimDemo2 , args=(dis_Gen,gLLsMat,notIn,gen1Name,gen2Name,0,300,gen_Gen,genList_GMat,conbamGen,dis_GenIndex,))

	p2 = multiprocessing.Process(target=getDisFunSimDemo2 , args=(dis_Gen,gLLsMat,notIn,gen1Name,gen2Name,300,lenth,gen_Gen,genList_GMat,conbamGen,dis_GenIndex,))

	# 运行多进程， 执行任务

	p1.start()

	p2.start()

	# 等待所有的子进程执行结束， 再执行主进程的内容

	p1.join()

	p2.join()

	print("任务执行结束.......")
	# disFuns = getDisFunSim(dis_Gen, gLLsMat, notIn, gen1Name, gen2Name)
	# np.savetxt(r'..\..\hmddV3\disfunS.txt', disFuns, fmt='%s', delimiter='\t', encoding='utf-8')

	print(time.process_time() - startTime)
