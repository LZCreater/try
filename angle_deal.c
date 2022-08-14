float GYROSCOPE_ANGLE_RATIO_X=1.03;              //X轴陀螺仪比例因子                    *****比例因子要改******
float GYROSCOPE_ANGLE_RATIO_Y =0.95;                 //Y轴陀螺仪比例因子
float GYROSCOPE_ANGLE_RATIO_Z=1;                 //Z轴陀螺仪比例因子
float AnglespeedOffset_X;                        //X轴零飘
float AnglespeedOffset_Y;                        //Y轴零飘
float AnglespeedOffset_Z;                        //Z轴零飘 **********（记得是浮点型）******
void get_anglespeed_bias()                      //获取三轴角速度零飘
{
  unsigned char i;
  for(i=0;i<100;i++)
  {
    Get_Gyro();    //获取一次陀螺仪最近储存的角速度
    systick_delay_ms(2);                            
    AnglespeedOffset_X+=mpu_gyro_x;
    AnglespeedOffset_Y+=mpu_gyro_y;
    AnglespeedOffset_Z+=mpu_gyro_z;
  }
    AnglespeedOffset_X /= 100;
    AnglespeedOffset_X /= 16.4;

    AnglespeedOffset_Y /= 100;
    AnglespeedOffset_Y /= 16.4;

    AnglespeedOffset_Z /= 100;
    AnglespeedOffset_Z /= 16.4;                 

}

int16 tempa,tempb,tempc,max,min;                 //用于陀螺仪滤波
int16 pre_Lastgyro_Y=0,Lastgyro_Y=0;             //Y轴陀螺仪均值滤波变量
extern float gyro_x, gyro_y, gyro_z;             
void get_anglespeed_y()                          //中值滤波获取y轴角速度
{
  Get_Gyro();                                   //获得一次角速度
//  gyro_y = mpu_gyro_y;
  tempa = pre_Lastgyro_Y;                     //陀螺仪均值滤波，保证得到陀螺仪值不跳变
  tempb = Lastgyro_Y;
  tempc = mpu_gyro_y;
  max = tempa > tempb ? tempa:tempb;
  max = max > tempc ? max:tempc;           
  min = tempa < tempb ? tempa:tempb;
  min = min < tempc ? min:tempc;          
  if(tempa > min  && tempa < max)   mpu_gyro_y = tempa;
  if(tempb > min  && tempb < max )  mpu_gyro_y = tempb;
  if(tempc > min  && tempc < max)   mpu_gyro_y = tempc; // 保证不跳变
  pre_Lastgyro_Y = Lastgyro_Y;//角速度递推赋值
  Lastgyro_Y = mpu_gyro_y;
//  mpu_gyro_y /= 16.4;    //这里再除以16.4 因为短整型int16是32767 然后除以2000的幅值        ((float)mpu_gyro_y)/16.4
  gyro_y = GYROSCOPE_ANGLE_RATIO_Y *(AnglespeedOffset_Y-((float)mpu_gyro_y)/16.4);//此处速度为正反馈，故零偏值减去AD值
  //陀螺仪（角速度）陀螺仪采集到的角速度归一化  
  //角速度归一化 得到的这个gyro_y才是最后真正得到的角速度 这里的比例因子 GYROSCOPE_ANGLE_RATIO_Y是重要 会影响角度融合相应的速度 
}

float k2=0.9,k1=0.05;

float x1,x2,y1;

float angle2,angle1;

float dt=0.04;                                 

float Erjielvbo(float angle_m,float gyro_m)     //二阶互补滤波 
{
    x1=(angle_m-angle2)*k2*k2;                  //0.9是常量（0.8）
    y1=y1+x1*dt;
    x2=y1+2*k2*(angle_m-angle2)+gyro_m;
    angle2=angle2+x2*dt;
    return angle2;
}

float Yijielvbo(float angle_m, float gyro_m)//采集后计算的角度和角加速度
{
  angle1=k1 * angle_m+ (1-k1) * (angle1 + gyro_m * dt);
  return angle1;
}

float angle, angle_dot;                         //角度和角速度

void Get_Angle(void)                            //获得当前角度
{
   get_anglespeed_y();                          //减零飘后的角速度
 //  Get_Gyro();
// gyro_y = mpu_gyro_y ;
   Get_AccData();                                
   Angle_A=atan2(mpu_acc_x,mpu_acc_z)*57.3;     
//   xg1 += gyro_y*dgt;           //累加一次 算比例因子       看加速度反正切求出的角度和角速度积分出来的角度是不是一样
   Kalman_Filter(Angle_A,gyro_y); 
   Angle1=angle;                                //卡尔曼滤波取最优角度  赋值
//   if(Angle1>90)Angle1=90;                      //限幅
//   if(Angle1<-90)Angle1=-90;
//   Angle2=Yijielvbo(Angle_A,gyro_y);            //一阶滤波
//   Angle3=Erjielvbo(Angle_A,gyro_y);            //二阶滤波  
}

 /*记载一下这里使用二阶互补滤波得到的角度变化反应很慢  */
/*
    下面是卡尔曼滤波 你直接复制整个就好
*/
float angle_0, angle_dot_0;//采集来的角度和角速度
float kdt=0.04;        //注意：dt的取值为kalman滤波器采样时间一下为运算中间变量
float P[2][2] = {{1,0},
                {0,1}};
float Pdot[4] ={0,0,0,0};
float Q_angle=0.001, Q_gyro=0.005; //角度数据置信度,角速度数据置信度
float R_angle=0.5 ,C_0 = 1; 
float q_bias, angle_err, PCt_0, PCt_1, E, K_0, K_1, t_0, t_1;
void Kalman_Filter(double angle_m,double gyro_m)
{
    angle+=(gyro_m-q_bias) * kdt;
    angle_err = angle_m - angle;
    Pdot[0]=Q_angle - P[0][1] - P[1][0];
    Pdot[1]=- P[1][1];
    Pdot[2]=- P[1][1];
    Pdot[3]=Q_gyro;
    P[0][0] += Pdot[0]*kdt;
    P[0][1] += Pdot[1]*kdt;
    P[1][0] += Pdot[2]*kdt;
    P[1][1] += Pdot[3]*kdt;
    PCt_0 = C_0 * P[0][0];
    PCt_1 = C_0 * P[1][0];
    E = R_angle + C_0 * PCt_0;
    K_0 = PCt_0 / E;
    K_1 = PCt_1 / E;
    t_0 = PCt_0;
    t_1 = C_0 * P[0][1];
    P[0][0] -= K_0 * t_0;
    P[0][1] -= K_0 * t_1;
    P[1][0] -= K_1 * t_0;
    P[1][1] -= K_1 * t_1;
    angle += K_0 * angle_err; //最优角度
    q_bias += K_1 * angle_err;
    angle_dot = gyro_m-q_bias;//最优角速度
}
