#include "HoleMesh.h"
#include <vector>
#include <queue>

#define gydebug

using namespace std;


int SingleHoleMesh::InitMesh(OMTriMesh & OriginMesh, int neighborRings, double distanceLimited) // exract
{
	holeMesh.clear();
	holeMesh.add_property(topDis);
	holeMesh.add_property(meshDis);
	queue<OpenMesh::VertexHandle> neighborVertices;
	OpenMesh::VertexHandle newV;
	OpenMesh::FaceHandle newFace;
	for (auto vit = VertBoundaryHs.begin(); vit != VertBoundaryHs.end(); vit++) {
		newV = holeMesh.add_vertex(OriginMesh.point(*vit));
		VertBoundaryHolemeshHs.push_back(newV);
		VertNew2OriId[newV.idx()] = vit->idx();
		VertOri2NewId[vit->idx()] = newV.idx();
		holeMesh.property(topDis, newV) = 0;
		holeMesh.property(meshDis, newV) = 0;
		neighborVertices.push(*vit);
	}
	int neibor=0;
	double distance,tempDis;
	OpenMesh::VertexHandle curV;
	while (neibor < neighborRings && !neighborVertices.empty())
	{
		curV = neighborVertices.front();
		neighborVertices.pop();

		newV = holeMesh.vertex_handle(VertOri2NewId[curV.idx()]);
		neibor = holeMesh.property(topDis, newV);
		distance = holeMesh.property(meshDis, newV);
		if (neibor == neighborRings - 1) {
			//need not this  neighbor
			break;
		}
		for (auto nVit = OriginMesh.vv_cwbegin(curV); nVit != OriginMesh.vv_cwend(curV); nVit++) {
			//if not in idx-map, 
				//if dis + newV.dis <distanceLimited and neibor< neighborRings
			//if in idx-map
				//if dis + newV.dis <distanceLimited 
				//if neibor< neighborRings
			if( VertOri2NewId.find( nVit->idx()) == VertOri2NewId.end()){
				tempDis = (OriginMesh.point(*nVit) - OriginMesh.point(curV)).length();
				if ( tempDis + distance < distanceLimited && neibor + 1 < neighborRings) {
					//add new point
					newV = holeMesh.add_vertex(OriginMesh.point(*nVit));
					VertNew2OriId[newV.idx()] = nVit->idx();
					VertOri2NewId[nVit->idx()] = newV.idx();
					holeMesh.property(topDis, newV) = neibor + 1;
					holeMesh.property(meshDis, newV) = tempDis + distance;
					neighborVertices.push(*nVit);
				}
			}
			else {
				tempDis = (OriginMesh.point(*nVit) - OriginMesh.point(curV)).length();
				if (tempDis + distance < holeMesh.property(meshDis, holeMesh.vertex_handle(VertOri2NewId[nVit->idx()]))) {
					//update new point
					holeMesh.property(meshDis, holeMesh.vertex_handle(VertOri2NewId[nVit->idx()])) = tempDis + distance;

				}
				if (neibor + 1 < holeMesh.property(topDis, holeMesh.vertex_handle(VertOri2NewId[nVit->idx()]))) {
					//update new point
					holeMesh.property(topDis, holeMesh.vertex_handle(VertOri2NewId[nVit->idx()])) = neibor + 1;
				}
			}

		}
		//add face
		//faces around curV
		for (auto nfit = OriginMesh.vf_cwbegin(curV); nfit != OriginMesh.vf_cwend(curV);nfit++) {
			if ( FacesOri2NewId.find(nfit->idx()) == FacesOri2NewId.end() ) {// has not been added into new mesh
				//check points of the face
				vector<OpenMesh::VertexHandle> triangleVertices;
				for (auto nFVit = OriginMesh.fv_cwbegin(*nfit); nFVit != OriginMesh.fv_cwend(*nfit); nFVit++) {
					if (VertOri2NewId.find(nFVit->idx()) == VertOri2NewId.end()) {//new vertex which is not in new mesh
						tempDis = (OriginMesh.point(*nFVit) - OriginMesh.point(curV)).length();
						//if (tempDis + distance < distanceLimited && neibor + 1 < neighborRings) {
						//add new point
						newV = holeMesh.add_vertex(OriginMesh.point(*nFVit));
						VertNew2OriId[newV.idx()] = nFVit->idx();
						VertOri2NewId[nFVit->idx()] = newV.idx();
						holeMesh.property(topDis, newV) = neibor + 1;
						holeMesh.property(meshDis, newV) = tempDis + distance;
						//neighborVertices.push(*nFVit);
						triangleVertices.push_back(newV);
					}
					else {
						triangleVertices.push_back( holeMesh.vertex_handle( VertOri2NewId.at( nFVit->idx() ) ) );
					}
				}
				if (triangleVertices.size() != 3) {
					std::cerr << "not a legal neighbor face" <<std::endl;
					continue;
				}
				newFace = holeMesh.add_face(triangleVertices[0], triangleVertices[2], triangleVertices[1]);//notice order
				FacesOri2NewId[nfit->idx()] = newFace.idx();
				FacesNew2OriId[newFace.idx()] = nfit->idx();
			}
		}
	}
	

#ifdef gydebugn
	string file = "G:\\gyGit\\data\\hole\\InitMesh"+to_string(AvergEdgeLen)+".obj";
	OpenMesh::IO::write_mesh(holeMesh, file);
	std::cout << "WriteMeshFile done" << file << std::endl;
#endif 
	MoveToCenter(CenterPoint);
#ifdef gydebugn
	file = "G:\\gyGit\\data\\hole\\InitMeshMove" + to_string(AvergEdgeLen) + ".obj";
	OpenMesh::IO::write_mesh(holeMesh, file);
	std::cout << "WriteMeshFile done" << file << std::endl;
#endif 
	return 0;
}

int SingleHoleMesh::Mergeback(OMTriMesh & OriginMesh, bool updateOldVertecis)
{
	MoveToCenter(-CenterPoint);
	OpenMesh::VertexHandle originV;
	OpenMesh::Vec3f	holeP;
	//vertece
	for  (OpenMesh::VertexHandle var : holeMesh.vertices() )
	{
		if ( VertNew2OriId.find(var.idx()) == VertNew2OriId.end()  ) {
			holeP = holeMesh.point(var);
			originV = OriginMesh.add_vertex(holeP);
			//add map
			VertNew2OriId[var.idx()] = originV.idx();
			VertOri2NewId[originV.idx()] = var.idx();
		}
		else if(updateOldVertecis){
			originV = OriginMesh.vertex_handle(VertNew2OriId[var.idx()]);
			holeP = holeMesh.point(var);
			OriginMesh.set_point(originV, holeP);
		}
		else {
			//do nothing
		}
	}
#ifdef gydebug
	string file = "G:\\gyGit\\data\\hole\\addV" + to_string(AvergEdgeLen) + ".obj";
	OpenMesh::IO::write_mesh(OriginMesh, file);
	std::cout << "WriteMeshFile done" << file << std::endl;
#endif 
	//face NewtFaceHs
	OpenMesh::VertexHandle originVs[3];
	OpenMesh::FaceHandle originF;
	for (OpenMesh::FaceHandle fa : holeMesh.faces())
	{
		if (FacesNew2OriId.find(fa.idx()) == FacesNew2OriId.end()) {
			int i = 0;
			for (auto v = holeMesh.fv_cwbegin(fa); v != holeMesh.fv_cwend(fa); v++,i++)
			{
				if (i >= 3) {
					i++;
					//return 3;//not a triangle
					break;
				}
				originVs[i] = OriginMesh.vertex_handle(VertNew2OriId[v->idx()]);
			}
			if (i > 3) {//after for it should be 3;after break it is 4
				//return 3;//not a triangle
				continue;
			}

			originF = OriginMesh.add_face(originVs[1], originVs[0], originVs[2]);
			//add map
			FacesNew2OriId[fa.idx()] = originF.idx();
			FacesOri2NewId[originF.idx()] = fa.idx();
		}
		else {
			//do nothing
		}
	}
#ifdef gydebug
	file = "G:\\gyGit\\data\\hole\\addF" + to_string(AvergEdgeLen) + ".obj";
	OpenMesh::IO::write_mesh(OriginMesh, file);
	std::cout << "WriteMeshFile done" << file << std::endl;
#endif 
	holeMesh.remove_property(topDis);
	holeMesh.remove_property(meshDis);
	return 0;
}

int SingleHoleMesh::FillHole(int type)
{
	if (0==type) {
		FillHoleCenterPoint();
	}
	else if (1==type) {
	}
	else {
		FillHoleCenterPoint();
	}
	

#ifdef gydebug
	string file = "G:\\gyGit\\data\\hole\\FillHole" + to_string(AvergEdgeLen) + ".obj";
	OpenMesh::IO::write_mesh(holeMesh, file);
	std::cout << "WriteMeshFile done" << file << std::endl;
#endif 
	return 0;
}

int SingleHoleMesh::FillHoleCenterPoint()
{
	//it is sure hole.VertHssize >= 3
	OpenMesh::VertexHandle CenterV = holeMesh.add_vertex(OpenMesh::Vec3f(0,0,0));//CenterPoint has been moved
	holeMesh.property(topDis, CenterV) = -1;

	//holeMesh.data(CenterV).type = -1;//new point
	auto startVertex = VertBoundaryHolemeshHs.begin();
	auto vIt = VertBoundaryHolemeshHs.begin();
	auto vItNext = vIt + 1;
	for (; vItNext != VertBoundaryHolemeshHs.end(); vIt++, vItNext++)
	{
		NewtFaceHs.push_back(holeMesh.add_face(*vIt, *vItNext, CenterV));//must be clock wise
	}
	NewtFaceHs.push_back(holeMesh.add_face(*vIt, *startVertex,  CenterV));
	return !(VertBoundaryHs.size() == NewtFaceHs.size());

}

int SingleHoleMesh::SmoothMesh(bool changeOrigin)
{
	//for each p
	//round by surrond

	//add to matrix

	//set value to point 

	return 0;
}

int SingleHoleMesh::MoveToCenter(OpenMesh::Vec3f center)
{
	//foreach p - CenterPoint
	for (auto v :holeMesh.vertices() ) {
		holeMesh.set_point(v, holeMesh.point(v) - center );
	}
	return 0;
}

///////////////////////////////////// HoleMeshHandle


int HoleMeshHandle::FindHoles()
{
	//OpenMesh::VPropHandleT<int> vertexType;
	//mesh.add_property(vertexType);
	holes.clear();
	for (auto vertexIt = mesh.vertices_begin(); vertexIt != mesh.vertices_end(); vertexIt++)
	{
		//find a boundary vertex, means a hole start here 
		if (mesh.is_boundary(*vertexIt) && mesh.data(*vertexIt).type < 1)//mesh.property(vertexType, *vertexIt) < 1)
		{
			//mesh.property(vertexType, *vertexIt) = 1;
			mesh.data(*vertexIt).type = 1;
			SingleHoleMesh aHole;
			aHole.AvergEdgeLen = 0;
			aHole.StartVertex = *vertexIt;
			aHole.VertBoundaryHs.push_back(*vertexIt);
			aHole.CenterPoint = mesh.point(*vertexIt);
			OpenMesh::VertexHandle nextVertex = *vertexIt;
			//along boundary push all boundary vertex
			do {
				//find a boundary neighbor vertex
				for (auto neighborIt = mesh.cvv_begin(nextVertex); neighborIt != mesh.cvv_end(nextVertex); neighborIt++)
				{
					if (*neighborIt == aHole.StartVertex)
					{
						nextVertex = *neighborIt; //break while
						break;
					}
					if (mesh.is_boundary(*neighborIt) && mesh.data(*neighborIt).type < 1)// mesh.property(vertexType, *neighborIt) < 1) //  wrong: vertex with more than one hole
					{
						//mesh.property(vertexType, *neighborIt) = 1;
						mesh.data(*neighborIt).type = 1;
						aHole.VertBoundaryHs.push_back(*neighborIt);
						aHole.CenterPoint += mesh.point(*neighborIt);
						aHole.AvergEdgeLen += (mesh.point(nextVertex) - mesh.point(*neighborIt)).length();
						nextVertex = *neighborIt;
						break;// assumpt only one  boundary vertex, so got find next
					}
					//else if (0 == mesh.data(*neighborIt).type) //neighbor ring 2
					//{
					//	mesh.data(*neighborIt).type = 2; // push to a queue;
					//	//debug
					//	mesh.set_point(*neighborIt, mesh.point(*neighborIt) + mesh.normal(*neighborIt) );
					//}
				}
				//
			} while (nextVertex != aHole.StartVertex && nextVertex.is_valid());

			if (aHole.VertBoundaryHs.size() < 3) {
				continue;//false hole
			}
			aHole.AvergEdgeLen /= (aHole.VertBoundaryHs.size() - 1);
			aHole.CenterPoint /= aHole.VertBoundaryHs.size();

			holes.push_back(aHole);
		}
	}
	return 0;
}


int HoleMeshHandle::ExtractHolesRegin(int neighborRings, double distanceLimited)
{
	//maybe do it in find
	FindHoles();
	//cout << "find holes:"<<holes.size() << endl;
	
	double farPointDis;
	//parallel
	for (vector<SingleHoleMesh>::iterator holeIt = holes.begin(); holeIt != holes.end();holeIt++) {
		farPointDis = holeIt->AvergEdgeLen * distanceLimited;
		holeIt->InitMesh(mesh, neighborRings+1, farPointDis);
	}
	return 0;
}

int HoleMeshHandle::FillHoles(int method)
{
	//parallel
	for (vector<SingleHoleMesh>::iterator holeIt = holes.begin(); holeIt != holes.end(); holeIt++) {
		holeIt->FillHole(method);
	}
	return 0;
}

int HoleMeshHandle::MergerBack(bool updateOldVertecis)
{
	//parallel
	for (vector<SingleHoleMesh>::iterator holeIt = holes.begin(); holeIt != holes.end(); holeIt++) {
		holeIt->Mergeback(mesh, updateOldVertecis);

#ifdef gydebug_merger
		string file = "G:\\gyGit\\data\\hole\\InitMesh.obj";
		OpenMesh::IO::write_mesh(holeMesh, file);
		std::cout << "WriteMeshFile done" << file << std::endl;
#endif 
	}

	return 0;
}
