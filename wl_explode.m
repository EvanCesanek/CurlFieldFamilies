classdef wl_explode < handle
%WL_EXPLODE is class created to simulate an exploding sphere
%
%   E = WL_EXPLODE(SIZE, ATTR, FRAG_SIZE) creates the explosion object with
%   radius of length SIZE, radius explosion fragments as FRAG_SIZE and
%   color of the fragments ATTR in [R, G, B] format.
%
%   E.ExplodePop([X, Y, Z]) starts explosion at [X, Y, Z] position 
%
%   E.ExplodeProcess(OBJ) processes the FSM for the explosion until the
%   duration has expired, OBJ is a pointer to the Experiment which is
%   passed to ExplodeProcess so that it can access sounds and files
%   required from the cfg property in OBJ.

    properties (Constant)
        EXPLODE_NONE  =   0;
        EXPLODE_WAITING = 1;
        EXPLODE_POP   =   2;
        EXPLODE_POPPING = 3;  
        EXPLODE_POPPED  = 4;
    end
    
    properties
        ExplodePosn;
        ExplodeSize;
        ExplodeFragments=16;
        ExplodeDuration=0.3;
        ExplodeSpeed=2.5;
        ExplodeFragmentSize;
        ExplodeAttr;
        ExplodeState=wl_explode.EXPLODE_NONE;
        ExplodeTimer;
        
    end
    
    methods
        function this=wl_explode(size, attr, frag_size, varargin)
            this.ExplodeSize = size;
            this.ExplodeAttr = attr;
            this.ExplodeFragmentSize = frag_size;
            this.ExplodeTimer = wl_timer;
            
            % optional initial properties 
            % variableSettings = ['ExplodeFragments' , 'ExplodeDuration', 'ExplodeSpeed', 'ExplodeState'];
            % Optional properties setter
            % for i=1:length(varargin)
            %     get(this, variableSettings(i)) = varargin{i};
            % end
            this.ExplodeState = this.EXPLODE_WAITING;
        end
        
        function ExplodeOpen(this, fragments, duration, speed)
            this.ExplodeFragments = fragments;
            this.ExplodeDuration = duration;
            this.ExplodeSpeed = speed;
            
        end
        
        function ExplodeClear(this)
            this.ExplodeState = this.EXPLODE_NONE;
        end
        
        function ExplodePop(this, pos)
            this.ExplodePosn = pos;
            if this.ExplodeState == this.EXPLODE_WAITING
                this.ExplodeState = this.EXPLODE_POP;
            end
        end
        
        function bool=ExplodePopping(this)
            bool = this.ExplodeState >= this.EXPLODE_POP && this.ExplodeState <= this.EXPLODE_POPPING;
        end
        
        function ExplodeProcess(this, obj)

            switch this.ExplodeState
                case this.EXPLODE_NONE
                    
                case this.EXPLODE_WAITING
                    
                case this.EXPLODE_POP
                    %wl_play_sound(obj, obj.cfg.explosion,obj.cfg.vol);
                    this.ExplodeTimer.Reset();
                    this.ExplodeState = this.ExplodeState + 1;
                    
                case this.EXPLODE_POPPING
                    
                    t = this.ExplodeTimer.GetTime;
                    if t >= this.ExplodeDuration
                        this.ExplodeState = this.ExplodeState + 1;
                    else
                        for i=0:this.ExplodeFragments-1
  
                            size = ((this.ExplodeDuration-t) / this.ExplodeDuration) * this.ExplodeFragmentSize;
                           
                            pos = this.ExplodePosn;
                            a = 2 * pi * i / this.ExplodeFragments;
                            
                            pos(1) = pos(1) + sin(a) * ((this.ExplodeSpeed*t) + (0.8*this.ExplodeSize));
                            pos(2) = pos(2) + cos(a) * ((this.ExplodeSpeed*t) + (0.8*this.ExplodeSize));

                            wl_draw_sphere(pos,size, this.ExplodeAttr,'Alpha',0.7);
                        end
                    end
                    
                case this.EXPLODE_POPPED
                    this.ExplodeState = this.EXPLODE_WAITING;
             
            end
        end
            
    end
    
end

