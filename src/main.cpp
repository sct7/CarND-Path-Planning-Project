#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "Eigen-3.3/Eigen/LU"
#include "json.hpp"
#include "spline.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2((map_y-y),(map_x-x));

	double angle = fabs(theta-heading);
  angle = min(2*pi() - angle, angle);

  if(angle > pi()/4)
  {
    closestWaypoint++;
    if (closestWaypoint == maps_x.size())
    {
      closestWaypoint = 0;
    }
  }

  return closestWaypoint;
}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}

vector<vector<double>> process_cars(const vector<vector<double>> &sensor_fusion, double car_s){
  //cout<<"2"<<endl;
  vector<vector<vector<double>>> sorted_cars;
  for (int i=0; i<3; i++) sorted_cars.push_back(vector<vector<double>>());

  
  //cout<<"3"<<endl;

  for (int c=0; c<sensor_fusion.size(); c++){
    //cout<<"car "<<c<<endl;
    vector<double> car = sensor_fusion[c];
    int lane = floor(car[6]/4);

    //cout<<lane<<endl;
    if (lane>=0 && lane<3){
      sorted_cars[lane].push_back(vector<double>());
      //cout<<sorted_cars[0].size()<<" "<<sorted_cars[1].size()<<" "<<sorted_cars[2].size()<<endl;
      for (int i=0; i<car.size(); i++){
        sorted_cars[lane][sorted_cars[lane].size()-1].push_back(car[i]);
      }
    }
    //cout<<"6"<<endl;
  }

  //cout<<"7"<<endl;

  vector<vector<double>> closest_cars;
  //for (int i=0; i<3; i++) sorted_cars.push_back(vector<double>());
  
  for (int lane = 0; lane<3; lane++){
    double closest_s = 10000;
    double closest_speed = 10000;
    for (int c = 0; c<sorted_cars[lane].size(); c++){
      vector<double> car = sorted_cars[lane][c];
      double check_s=car[5];
      double dist = fmod(check_s-car_s+6945.554, 6945.554);
      if (dist<closest_s){
        closest_s=dist;
        double check_xv = car[3];
        double check_yv = car[4];
        double check_speed = sqrt(pow(check_xv,2)+pow(check_yv,2))*2.24;
        closest_speed = check_speed;
      }
    }
    //cout<<closest_s<<" "<<closest_speed<<endl;
    vector<double> c_car = {closest_s, closest_speed};
    closest_cars.push_back(c_car);
  }
  return closest_cars;
}

bool is_clear(const vector<vector<double>> &sensor_fusion, double end_path_s, int target_lane){
  for (int c=0; c<sensor_fusion.size(); c++){
    vector<double> car = sensor_fusion[c];
    int check_lane = floor(car[6]/4);
    //cout<<check_lane<<endl;
    if (check_lane==target_lane){

      double check_s=car[5];
      double check_xv = car[3];
      double check_yv = car[4];
      double check_speed = sqrt(pow(check_xv,2)+pow(check_yv,2)); // m/s NOT MPH
      //cout<<"checking a car here "<<end_path_s<<" "<<check_s<<" "<<check_speed<<" "<<fabs(check_s+check_speed-end_path_s)<<endl;
      
      if (fabs(check_s+check_speed-end_path_s)<7){
        return false;
      }
    }
  }
  return true;
}

int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }

  double v_ref = 0.224;

  h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy, &v_ref](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") {
          // j[1] is the data JSON object
          
        	// Main car's localization Data
          	double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
          	double car_yaw = j[1]["yaw"];
          	double car_speed = j[1]["speed"];

          	// Previous path data given to the Planner
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
          	// Previous path's end s and d values 
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	auto sensor_fusion = j[1]["sensor_fusion"];

            int prev_size = previous_path_x.size();

          	json msgJson;

            int lane = floor(end_path_d/4);
            
            double TARGET_V =49.5;
            int FUTURE_POINTS = 50;


            if (prev_size>0){

              vector<vector<double>> cars = sensor_fusion;
              //cout<<"1"<<endl;
              vector<vector<double>> closest_cars = process_cars(cars, car_s);
              //cout<<"9"<<endl;

              

              int goal_lane=lane;
              int target_lane = lane;
              

              double lane_dist = closest_cars[lane][0];
              double lane_speed = closest_cars[lane][1];

              
                //find the fastest lane
                double max_speed=0;
                for (int l=0; l<3; l++){
                  double candidate_speed = min(closest_cars[l][1], TARGET_V);
                  if (closest_cars[l][0]>30) candidate_speed = TARGET_V + int(l==lane); //break ties by sticking w/ lane
                  if (candidate_speed>max_speed){
                      max_speed = candidate_speed;
                      goal_lane=l;
                  }
                }

              if (goal_lane>lane){
                cout<<"Attempting right ";
                if (is_clear(sensor_fusion, end_path_s, lane+1)){
                  cout<<"Right clear"<<endl;
                  target_lane+=1;
                }else{
                  cout<<"Right NOT clear"<<endl;
                }
              }else if (goal_lane<lane){
                cout<<"Attempting left ";
                if (is_clear(sensor_fusion, end_path_s, lane-1)){
                  cout<<"Left clear"<<endl;
                  target_lane-=1;
                }else{
                  cout<<"Left NOT clear"<<endl;
                }
              }
              cout<<goal_lane<<" "<<closest_cars[0][0]<<" "<<closest_cars[0][1]<<" "<<closest_cars[1][0]<<" "<<closest_cars[1][1]<<" "<<closest_cars[2][0]<<" "<<closest_cars[2][1]<<endl;
              //cout<<target_lane<<" "<<lane<<" "<<lane_speed<<endl;

              if (lane_dist<14 && lane_speed<v_ref){
                cout<<"slowing down"<<endl;
                v_ref-=0.224;
              }else if (v_ref<TARGET_V){
                v_ref+=0.224;
              }

              lane=target_lane;
            }else{
              end_path_s=car_s;
            }


            vector<double> x_pts;
            vector<double> y_pts;

            double ref_x=car_x;
            double ref_y=car_y;
            double ref_yaw = car_yaw;

            

            if (prev_size<2){
              //cout<<"sup"<<endl;
              x_pts.push_back(car_x-cos(car_yaw));
              y_pts.push_back(car_y-sin(car_yaw));
              x_pts.push_back(car_x);
              y_pts.push_back(car_y);
              
            }else{
              //cout<<"sup 2"<<endl;
              ref_x = previous_path_x[prev_size-1];
              ref_y = previous_path_y[prev_size-1];

              double prev_2_x = previous_path_x[prev_size-2];
              double prev_2_y = previous_path_y[prev_size-2];

              ref_yaw = atan2((ref_y-prev_2_y),(ref_x-prev_2_x));

              x_pts.push_back(prev_2_x);
              y_pts.push_back(prev_2_y);
              x_pts.push_back(ref_x);
              y_pts.push_back(ref_y);
            }



            for (int i=1; i<4; i++){
              auto next_xy = getXY(end_path_s+30*i, 2+lane*4, map_waypoints_s, map_waypoints_x, map_waypoints_y);
              x_pts.push_back(next_xy[0]);
              y_pts.push_back(next_xy[1]);
            }

            //cout<<"$$$$$$"<<endl;
            //ROTATE EVERYTHING INTO CAR FRAME
            for (int i=0; i<x_pts.size(); i++){
             // cout<<"--"<<endl;
              //cout<<x_pts[i]<<" "<<y_pts[i]<<endl;
              double shift_x = x_pts[i]-ref_x;
              double shift_y = y_pts[i]-ref_y;

              x_pts[i] = shift_x*cos(0-ref_yaw)-shift_y*sin(0-ref_yaw);
              y_pts[i] = shift_y*cos(0-ref_yaw)+shift_x*sin(0-ref_yaw);

              //cout<<x_pts[i]<<" "<<y_pts[i]<<endl;
            }

            tk::spline s;
            
            s.set_points(x_pts, y_pts);

            vector<double> next_x_vals;
            vector<double> next_y_vals;

            //THESE ARE IN GLOBAL FRAME
            for (int i=0; i<previous_path_x.size(); i++){
              next_x_vals.push_back(previous_path_x[i]);
              next_y_vals.push_back(previous_path_y[i]);
            }

            double t_x = 30;
            double t_y = s(t_x);
            double dist = distance(0, 0, t_x, t_y);
            double N = dist/(v_ref/2.24*0.02);

            double x_gap = t_x/N;
            //cout<<t_x<<" "<<t_y<<" "<<dist<<" "<<N<<" "<<v_ref<<endl;
            //cout<<"x_gap "<<x_gap<<endl;

            //cout<<"adding "<<FUTURE_POINTS-prev_size<<" points"<<endl;
            for (int i=1; i<=(FUTURE_POINTS-prev_size); i++){

              double new_x = x_gap*i;
              double new_y = s(x_gap*i);

              //cout<<new_x<<" "<<new_y<<endl;

              //ROTATE BACK INTO GLOBAL
              double temp_x = new_x;
              double temp_y = new_y;

              new_x = temp_x*cos(ref_yaw)-temp_y*sin(ref_yaw);
              new_y = temp_y*cos(ref_yaw)+temp_x*sin(ref_yaw);

              new_x+=ref_x;
              new_y+=ref_y;

              next_x_vals.push_back(new_x);
              next_y_vals.push_back(new_y);
            }

            
            //cout<<"all 50 points"<<endl;
            for (int i=0; i<next_x_vals.size(); i++){
              //cout<<next_x_vals[i]<<" "<<next_y_vals[i]<<endl;
            }
            //cout<<"end"<<endl;
            

          	msgJson["next_x"] = next_x_vals;
          	msgJson["next_y"] = next_y_vals;

          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconds(1000));
          	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
          
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
