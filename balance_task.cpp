/*============================================================================
==============================================================================
                      
                              balance_task.c
 
==============================================================================
Remarks:

      sekeleton to create the sample task

============================================================================*/

// system headers
#include "SL_system_headers.h"

/* SL includes */
#include "SL.h"
#include "SL_user.h"
#include "SL_tasks.h"
#include "SL_task_servo.h"
#include "SL_kinematics.h"
#include "SL_dynamics.h"
#include "SL_collect_data.h"
#include "SL_shared_memory.h"
#include "SL_man.h"

// defines

// local variables
static double      start_time = 0.0;
static SL_DJstate  target[N_DOFS+1];
static SL_Cstate   cog_target;
static SL_Cstate   cog_traj;
static SL_Cstate   cog_ref;
static double      delta_t = 0.01;
static double      duration = 5.0;//10.0; //used to be 10
static double      time_to_go;
static int         which_step;
// possible states of a state machine
enum Steps {
  ASSIGN_COG_TARGET_R,
  MOVE_TO_COG_TARGET_R,
  ASSIGN_JOINT_TARGET_LIFT_UP_R,
  MOVE_JOINT_TARGET_LIFT_UP_R,
  ASSIGN_JOINT_TARGET_LIFT_DOWN_R,
  MOVE_JOINT_TARGET_LIFT_DOWN_R,
  ASSIGN_COG_TARGET_L,
  MOVE_TO_COG_TARGET_L,
  ASSIGN_JOINT_TARGET_LIFT_UP_L,
  MOVE_JOINT_TARGET_LIFT_UP_L,
  ASSIGN_JOINT_TARGET_LIFT_DOWN_L,
  MOVE_JOINT_TARGET_LIFT_DOWN_L
};

// variables for COG control
static iMatrix     stat;
static Matrix      Jccogp;
static Matrix      NJccog;
static Matrix      fc;

// global functions 
extern "C" void
add_balance_task( void );

// local functions
static int  init_balance_task(void);
static int  run_balance_task(void);
static int  change_balance_task(void);

static int 
min_jerk_next_step (double x,double xd, double xdd, double t, double td, double tdd,
		    double t_togo, double dt,
		    double *x_next, double *xd_next, double *xdd_next);


/*****************************************************************************
******************************************************************************
Function Name	: add_balance_task
Date		: Feb 1999
Remarks:

adds the task to the task menu

******************************************************************************
Paramters:  (i/o = input/output)

none

*****************************************************************************/
void
add_balance_task( void )
{
  int i, j;
  
  addTask("Balance Task", init_balance_task, 
	  run_balance_task, change_balance_task);

}    

/*****************************************************************************
******************************************************************************
  Function Name	: init_balance_task
  Date		: Dec. 1997

  Remarks:

  initialization for task

******************************************************************************
  Paramters:  (i/o = input/output)

       none

 *****************************************************************************/
static int 
init_balance_task(void)
{
  int j, i;
  int ans;
  static int firsttime = TRUE;
  
  if (firsttime){
    firsttime = FALSE;

    // allocate memory
    stat   = my_imatrix(1,N_ENDEFFS,1,2*N_CART);
    Jccogp = my_matrix(1,N_DOFS,1,N_CART);
    NJccog = my_matrix(1,N_DOFS,1,N_DOFS+2*N_CART);
    fc     = my_matrix(1,N_ENDEFFS,1,2*N_CART);
    // this is an indicator which Cartesian components of the endeffectors are constraints
    // i.e., both feet are on the ground and cannot move in position or orientation
    stat[RIGHT_FOOT][1] = TRUE;
    stat[RIGHT_FOOT][2] = TRUE;
    stat[RIGHT_FOOT][3] = TRUE;
    stat[RIGHT_FOOT][4] = TRUE;
    stat[RIGHT_FOOT][5] = TRUE;
    stat[RIGHT_FOOT][6] = TRUE;

    stat[LEFT_FOOT][1] = TRUE;
    stat[LEFT_FOOT][2] = TRUE;
    stat[LEFT_FOOT][3] = TRUE;
    stat[LEFT_FOOT][4] = TRUE;
    stat[LEFT_FOOT][5] = TRUE;
    stat[LEFT_FOOT][6] = TRUE;

  }

  // prepare going to the default posture
  bzero((char *)&(target[1]),N_DOFS*sizeof(target[1]));
  for (i=1; i<=N_DOFS; i++)
    target[i] = joint_default_state[i];

  // go to the target using inverse dynamics (ID)
  if (!go_target_wait_ID(target)) 
    return FALSE;

  // ready to go
  ans = 999;
  while (ans == 999) {
    if (!get_int("Enter 1 to start or anthing else to abort ...",ans,&ans))
      return FALSE;
  }
  
  // only go when user really types the right thing
  if (ans != 1) 
    return FALSE;

  start_time = task_servo_time;
  printf("start time = %.3f, task_servo_time = %.3f\n", 
	 start_time, task_servo_time);

  // start data collection
  scd();

  // state machine starts at ASSIGN_COG_TARGET
  which_step = ASSIGN_COG_TARGET_R;

  return TRUE;
}

/*****************************************************************************
******************************************************************************
  Function Name	: run_balance_task
  Date		: Dec. 1997

  Remarks:

  run the task from the task servo: REAL TIME requirements!

******************************************************************************
  Paramters:  (i/o = input/output)

  none

 *****************************************************************************/
static int 
run_balance_task(void)
{
  int j, i, n;

  double task_time;
  double kp = 0.1;
  // ******************************************
  // NOTE: all array indices start with 1 in SL
  // ******************************************

  task_time = task_servo_time - start_time;

  // the following code computes the contraint COG Jacobian 
  // Jccogp is an N_DOFS x N_CART matrix
  // NJccog is an N_DOFS x N_DOF+2*N_CART matrix -- most likely this is not needed

  compute_cog_kinematics(stat, TRUE, FALSE, TRUE, Jccogp, NJccog);
  // switch according to the current state of the state machine
  switch (which_step) {
    
  case ASSIGN_COG_TARGET_R:

    // what is the target for the COG?
    bzero((void *)&cog_target,sizeof(cog_target));
    cog_target.x[_X_] = 0.055; // pos: right
    cog_target.x[_Y_] = 0.015; // pos: forward
    cog_target.x[_Z_] = -.115;

    // the structure cog_des has the current position of the COG computed from the
    // joint_des_state of the robot. cog_des should track cog_traj
    bzero((void *)&cog_traj,sizeof(cog_traj));
    for (i=1; i<=N_CART; ++i)
      cog_traj.x[i] = cog_des.x[i];

    // time to go
    time_to_go = duration;

    // switch to next step of state machine
    which_step = MOVE_TO_COG_TARGET_R;
    //printf("COG assigned well\n");
    // Question: WHY change the state right away? How to make sure time_to_go = 0?
    // Ans: no need to wait after initialization
    break;

  case MOVE_TO_COG_TARGET_R: // this is for inverse kinematics control

    // plan the next step of cog with min jerk
    for (i=1; i<=N_CART; ++i) {
      min_jerk_next_step(cog_traj.x[i],
			 cog_traj.xd[i],
			 cog_traj.xdd[i],
			 cog_target.x[i],
			 cog_target.xd[i],
			 cog_target.xdd[i],
			 time_to_go,
			 delta_t,
			 &(cog_traj.x[i]),
			 &(cog_traj.xd[i]),
			 &(cog_traj.xdd[i]));
    }

    // inverse kinematics: we use a P controller to correct for tracking erros
    for (i=1; i<=N_CART; ++i)
      cog_ref.xd[i] = kp*(cog_traj.x[i] - cog_des.x[i]) + cog_traj.xd[i]; // why only xd?
/*    
  if (!inverseKinematics(target,endeff,joint_opt_state,&cog_ref,stat,delta_t)) {
//    printf("IK Condition number = 0\n"); // debug
    freeze();
    return FALSE;
  }
*/
    // compute the joint_des_state[i].th and joint_des_state[i].thd  
    for (i=1; i<=N_DOFS; ++i) {
      // intialize to zero
      joint_des_state[i].thd  = 0;
      joint_des_state[i].thdd = 0;
      joint_des_state[i].uff  = 0;

// dw 4/8 modified
      for (j=1; j<=N_CART; ++j) {
						  joint_des_state[i].th += delta_t * Jccogp[i][j] * cog_ref.xd[j];
        joint_des_state[i].thd += Jccogp[i][j] * cog_ref.xd[j];
      }        
      //joint_des_state[i].thdd = target[i].thdd; 
/* didn't help
//added		
      if (joint_des_state[i].th > joint_range[i][MAX_THETA]) {
        joint_des_state[i].th = joint_range[i][MAX_THETA];
//        joint_des_state[i].thd = 0.0;
//        joint_des_state[i].thdd = 0.0;
      }
      if (joint_des_state[i].th < joint_range[i][MIN_THETA]) {
        joint_des_state[i].th = joint_range[i][MIN_THETA];
//        joint_des_state[i].thd = 0.0;
//        joint_des_state[i].thdd = 0.0;
      }
*/						
    }

    // decrement time to go
    time_to_go -= delta_t;
    //printf("COG moved well\n");
    if (time_to_go <= 0) {
        which_step = ASSIGN_JOINT_TARGET_LIFT_UP_R;
        time_to_go = duration;
    }

    break;

  case ASSIGN_JOINT_TARGET_LIFT_UP_R:

    // initialize the target structure from the joint_des_state
    for (i=1; i<=N_DOFS; ++i)
      target[i] = joint_des_state[i];

    target[L_HFE].th = 0.6;
    target[L_KFE].th = 1.2;
    target[L_AFE].th = 0.6;// neg means go down

    time_to_go = duration;  // this may be too fast -- maybe a slower movement is better

    which_step = MOVE_JOINT_TARGET_LIFT_UP_R;
    //printf("lift left leg up assigned well\n");
    break;

  case MOVE_JOINT_TARGET_LIFT_UP_R:

    // compute the update for the desired states
    for (i=1; i<=N_DOFS; ++i) {
      min_jerk_next_step(joint_des_state[i].th,
			 joint_des_state[i].thd,
			 joint_des_state[i].thdd,
			 target[i].th,
			 target[i].thd,
			 target[i].thdd,
			 time_to_go,
			 delta_t,
			 &(joint_des_state[i].th),
			 &(joint_des_state[i].thd),
			 &(joint_des_state[i].thdd));
    }

    // decrement time to go
    time_to_go -= delta_t;
//    printf("In MOVE_JOINT_TARGET_LIFT_UP\n");
    if (time_to_go <= 0){
        which_step = ASSIGN_JOINT_TARGET_LIFT_DOWN_R;
        time_to_go = duration;
    }
    //printf("lift left leg up moved well\n");
    break;

  case ASSIGN_JOINT_TARGET_LIFT_DOWN_R:
      
      // initialize the target structure from the joint_des_state
      for (i=1; i<=N_DOFS; ++i)
          target[i] = joint_des_state[i];
      
      target[L_HFE].th = 0.05;
      target[L_KFE].th = 0.002;
      target[L_AFE].th = -0.04;// neg means go down
      
      time_to_go = duration;  // this may be too fast -- maybe a slower movement is better
      
      which_step = MOVE_JOINT_TARGET_LIFT_DOWN_R;
      
      //printf("lift left leg downb assigned well\n");
      break;
      
  case MOVE_JOINT_TARGET_LIFT_DOWN_R:
      
      // compute the update for the desired states
      for (i=1; i<=N_DOFS; ++i) {
          min_jerk_next_step(joint_des_state[i].th,
                             joint_des_state[i].thd,
                             joint_des_state[i].thdd,
                             target[i].th,
                             target[i].thd,
                             target[i].thdd,
                             time_to_go,
                             delta_t,
                             &(joint_des_state[i].th),
                             &(joint_des_state[i].thd),
                             &(joint_des_state[i].thdd));
      }
      
      // decrement time to go
      time_to_go -= delta_t;
      //    printf("In MOVE_JOINT_TARGET_LIFT_UP\n");
      if (time_to_go <= 0){
          which_step = ASSIGN_COG_TARGET_L;
          time_to_go = duration;
      }
     //printf("lift left leg down moved well\n");
    break;
    
  case ASSIGN_COG_TARGET_L:
      
      // what is the target for the COG?
      bzero((void *)&cog_target,sizeof(cog_target));
      cog_target.x[_X_] = -0.055; // pos: right
      cog_target.x[_Y_] = 0.015; // pos: forward
      cog_target.x[_Z_] = .115;
      
      // the structure cog_des has the current position of the COG computed from the
      // joint_des_state of the robot. cog_des should track cog_traj
      bzero((void *)&cog_traj,sizeof(cog_traj));
      for (i=1; i<=N_CART; ++i)
          cog_traj.x[i] = cog_des.x[i];
      
      // time to go
      time_to_go = duration;
      
      // switch to next step of state machine
      which_step = MOVE_TO_COG_TARGET_L;
      
      // Question: WHY change the state right away? How to make sure time_to_go = 0?
      // Ans: no need to wait after initialization
      //printf("COG assigned well\n");
      break;
      
  case MOVE_TO_COG_TARGET_L: // this is for inverse kinematics control
      
      // plan the next step of cog with min jerk
      for (i=1; i<=N_CART; ++i) {
          min_jerk_next_step(cog_traj.x[i],
                             cog_traj.xd[i],
                             cog_traj.xdd[i],
                             cog_target.x[i],
                             cog_target.xd[i],
                             cog_target.xdd[i],
                             time_to_go,
                             delta_t,
                             &(cog_traj.x[i]),
                             &(cog_traj.xd[i]),
                             &(cog_traj.xdd[i]));
      }
      
      // inverse kinematics: we use a P controller to correct for tracking erros
      for (i=1; i<=N_CART; ++i)
          cog_ref.xd[i] = kp*(cog_traj.x[i] - cog_des.x[i]) + cog_traj.xd[i]; // why only xd?
      /*
       if (!inverseKinematics(target,endeff,joint_opt_state,&cog_ref,stat,delta_t)) {
       //    printf("IK Condition number = 0\n"); // debug
       freeze();
       return FALSE;
       }
       */
      // compute the joint_des_state[i].th and joint_des_state[i].thd
      for (i=1; i<=N_DOFS; ++i) {
          // intialize to zero
          joint_des_state[i].thd  = 0;
          joint_des_state[i].thdd = 0;
          joint_des_state[i].uff  = 0;
          
          // dw 4/8 modified
          for (j=1; j<=N_CART; ++j) {
              joint_des_state[i].th += delta_t * Jccogp[i][j] * cog_ref.xd[j];
              joint_des_state[i].thd += Jccogp[i][j] * cog_ref.xd[j];
          }
          //joint_des_state[i].thdd = target[i].thdd;
          /* didn't help
           //added
           if (joint_des_state[i].th > joint_range[i][MAX_THETA]) {
           joint_des_state[i].th = joint_range[i][MAX_THETA];
           //        joint_des_state[i].thd = 0.0;
           //        joint_des_state[i].thdd = 0.0;
           }
           if (joint_des_state[i].th < joint_range[i][MIN_THETA]) {
           joint_des_state[i].th = joint_range[i][MIN_THETA];
           //        joint_des_state[i].thd = 0.0;
           //        joint_des_state[i].thdd = 0.0;
           }
           */
      }
      
      // decrement time to go
      time_to_go -= delta_t;
      if (time_to_go <= 0) {
          which_step = ASSIGN_JOINT_TARGET_LIFT_UP_L;
          time_to_go = duration;
          //freeze(); // added for now
      }
      //printf("COG moved well\n");
      break;
      
  case ASSIGN_JOINT_TARGET_LIFT_UP_L:
      
      // initialize the target structure from the joint_des_state
      for (i=1; i<=N_DOFS; ++i)
          target[i] = joint_des_state[i];
      
      target[R_HFE].th = 0.6;
      target[R_KFE].th = 1.2;
      target[R_AFE].th = 0.6;// neg means go down
      
      time_to_go = duration;  // this may be too fast -- maybe a slower movement is better
      
      which_step = MOVE_JOINT_TARGET_LIFT_UP_L;
      
      //    time_to_go -= delta_t;
      //    if (time_to_go <= 0) {
      //      which_step = MOVE_JOINT_TARGET_LIFT_UP;
      //                        time_to_go = duration;
      //                        //freeze(); // added for now
      //    }
   // printf("lift right leg up assigned well\n");
      break;
      
  case MOVE_JOINT_TARGET_LIFT_UP_L:
      
      // compute the update for the desired states
      for (i=1; i<=N_DOFS; ++i) {
          min_jerk_next_step(joint_des_state[i].th,
                             joint_des_state[i].thd,
                             joint_des_state[i].thdd,
                             target[i].th,
                             target[i].thd,
                             target[i].thdd,
                             time_to_go,
                             delta_t,
                             &(joint_des_state[i].th),
                             &(joint_des_state[i].thd),
                             &(joint_des_state[i].thdd));
      }
      
      // decrement time to go
      time_to_go -= delta_t;
      //    printf("In MOVE_JOINT_TARGET_LIFT_UP\n");
      if (time_to_go <= 0){
          which_step = ASSIGN_JOINT_TARGET_LIFT_DOWN_L;
          time_to_go = duration;
      }
    //printf("lift right leg up moved well\n");
      break;
      
  case ASSIGN_JOINT_TARGET_LIFT_DOWN_L:
      
      // initialize the target structure from the joint_des_state
      for (i=1; i<=N_DOFS; ++i)
          target[i] = joint_des_state[i];
      
      target[R_HFE].th = -0.05;
      target[R_KFE].th = 0.002;
      target[R_AFE].th = -0.04;// neg means go down
      
      time_to_go = duration;  // this may be too fast -- maybe a slower movement is better
      
      which_step = MOVE_JOINT_TARGET_LIFT_DOWN_L;
      
      //    time_to_go -= delta_t;
      //    if (time_to_go <= 0) {
      //      which_step = MOVE_JOINT_TARGET_LIFT_UP;
      //                        time_to_go = duration;
      //                        //freeze(); // added for now
      //    }
      
      break;
      
  case MOVE_JOINT_TARGET_LIFT_DOWN_L:
      
      // compute the update for the desired states
      for (i=1; i<=N_DOFS; ++i) {
          min_jerk_next_step(joint_des_state[i].th,
                             joint_des_state[i].thd,
                             joint_des_state[i].thdd,
                             target[i].th,
                             target[i].thd,
                             target[i].thdd,
                             time_to_go,
                             delta_t,
                             &(joint_des_state[i].th),
                             &(joint_des_state[i].thd),
                             &(joint_des_state[i].thdd));
      }
      
      // decrement time to go
      time_to_go -= delta_t;
      //    printf("In MOVE_JOINT_TARGET_LIFT_UP\n");
      if (time_to_go <= 0){
          //printf("Done\n");
          freeze();
      }
				
  }

  // this is a special inverse dynamics computation for a free standing robot
//  inverseDynamicsFloat(delta_t, stat, TRUE, joint_des_state, NULL, NULL, fc);


  return TRUE;
}

/*****************************************************************************
******************************************************************************
  Function Name	: change_balance_task
  Date		: Dec. 1997

  Remarks:

  changes the task parameters

******************************************************************************
  Paramters:  (i/o = input/output)

  none

 *****************************************************************************/
static int 
change_balance_task(void)
{
  int    ivar;
  double dvar;

  get_int("This is how to enter an integer variable",ivar,&ivar);
  get_double("This is how to enter a double variable",dvar,&dvar);

  return TRUE;

}


/*!*****************************************************************************
 *******************************************************************************
\note  min_jerk_next_step
\date  April 2014
   
\remarks 

Given the time to go, the current state is updated to the next state
using min jerk splines

 *******************************************************************************
 Function Parameters: [in]=input,[out]=output

 \param[in]          x,xd,xdd : the current state, vel, acceleration
 \param[in]          t,td,tdd : the target state, vel, acceleration
 \param[in]          t_togo   : time to go until target is reached
 \param[in]          dt       : time increment
 \param[in]          x_next,xd_next,xdd_next : the next state after dt

 ******************************************************************************/
static int 
min_jerk_next_step (double x,double xd, double xdd, double t, double td, double tdd,
		    double t_togo, double dt,
		    double *x_next, double *xd_next, double *xdd_next)

{
  double t1,t2,t3,t4,t5;
  double tau,tau1,tau2,tau3,tau4,tau5;
  int    i,j;

  // a safety check
  if (dt > t_togo || dt <= 0) {
    return FALSE;
  }

  t1 = dt;
  t2 = t1 * dt;
  t3 = t2 * dt;
  t4 = t3 * dt;
  t5 = t4 * dt;

  tau = tau1 = t_togo;
  tau2 = tau1 * tau;
  tau3 = tau2 * tau;
  tau4 = tau3 * tau;
  tau5 = tau4 * tau;

  // calculate the constants
  const double dist   = t - x;
  const double p1     = t;
  const double p0     = x;
  const double a1t2   = tdd;
  const double a0t2   = xdd;
  const double v1t1   = td;
  const double v0t1   = xd;
  
  const double c1 = 6.*dist/tau5 + (a1t2 - a0t2)/(2.*tau3) - 
    3.*(v0t1 + v1t1)/tau4;
  const double c2 = -15.*dist/tau4 + (3.*a0t2 - 2.*a1t2)/(2.*tau2) +
    (8.*v0t1 + 7.*v1t1)/tau3; 
  const double c3 = 10.*dist/tau3+ (a1t2 - 3.*a0t2)/(2.*tau) -
    (6.*v0t1 + 4.*v1t1)/tau2; 
  const double c4 = xdd/2.;
  const double c5 = xd;
  const double c6 = x;
  
  *x_next   = c1*t5 + c2*t4 + c3*t3 + c4*t2 + c5*t1 + c6;
  *xd_next  = 5.*c1*t4 + 4*c2*t3 + 3*c3*t2 + 2*c4*t1 + c5;
  *xdd_next = 20.*c1*t3 + 12.*c2*t2 + 6.*c3*t1 + 2.*c4;
  
  return TRUE;
}
