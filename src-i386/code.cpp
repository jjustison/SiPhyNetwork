#include <Rcpp.h>
using namespace Rcpp;

IntegerVector which2(const LogicalVector & x) {
  IntegerVector v = seq(0, x.size()-1);
  v= v[x];
  return(v);
}

int which1(const LogicalVector & x) {
  IntegerVector v = seq(0, x.size()-1);
  v= v[x];
  int w = v[0];
  return(w);
}

void biconnectedHelper(const IntegerMatrix &edges,const int nd,int &index, IntegerVector &edgeS, List &blobs, IntegerVector &depth, IntegerVector &lowpoint, IntegerVector &prev){
  int children = 0;

  depth[nd-1]    = index;
  lowpoint[nd-1] = index;
  index++;
  IntegerVector adj_edges = which2( (edges(_,0)==nd) | (edges(_,1)==nd) );
  for(auto e : adj_edges){
    IntegerVector nd_edge=edges(e,_);
    IntegerVector child_nd = nd_edge[nd_edge != nd];

    int w = child_nd[0];
    if(depth[w-1]==-1){
      prev[w-1]= nd;
      children++;
      edgeS.push_back(e);
      biconnectedHelper(edges,w,index,edgeS,blobs,depth,lowpoint,prev);

      lowpoint[nd-1]= min(as<IntegerVector>(lowpoint[nd_edge-1]));

      if(( (prev[nd-1]==-1) & (children >1) ) | ((prev[nd-1] != -1) & (lowpoint[w-1] >= depth[nd-1]) )) { //The node is an articulation point or the root
        //We make a new blob
        IntegerVector indices = seq(0, edgeS.size()-1);
        int start_ind = which1(edgeS == e); //The blob are all edges including and past this one
        if(start_ind != (indices.size()-1)){//only add the blob if it has more than one edge in it
			    blobs.push_back(as<IntegerVector>(edgeS[indices>=start_ind]));//Add the blob
        }
        edgeS = edgeS[indices<start_ind]; //Remove the blob edges from edgeS
      }
    } else if( (prev[nd-1]!=w) & (depth[nd-1]>depth[w-1]) ){ //e is a back edge
        lowpoint[nd-1] = std::min(lowpoint[nd-1],depth[w-1]);
        edgeS.push_back(e);
    }
  }
}

//' Find the biconnected components of a phylogeny
//'
//' @title Biconnected Components
//' @description Find the biconnected components of a phylogeny
//' @param edges The `edge` matrix of a `phylo` object.
//' @param rt The root node of the phylogeny
//' @param nNode The number of nodes in the tree
//' @details biconnected components
//' @return A list with containing a vector of nodes for each biconnected component
//' @export
// [[Rcpp::export]]
List biconnectedComponents(IntegerMatrix edges,int rt,int nNode){
  IntegerVector depth (nNode,-1);
  IntegerVector lowpoint (nNode,-1);
  IntegerVector prev(nNode,-1);

  IntegerVector edgeS;
  List blobs;
  int index = 0;

  biconnectedHelper(edges, rt, index,edgeS,blobs,depth,lowpoint,prev);

  if(edgeS.size()>1){
    blobs.push_back(edgeS);
  }
  return(blobs);
}


