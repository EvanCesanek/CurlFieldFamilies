function wl_error_state(WL, msg, state)
%WL_ERROR_STATE plays sounds corresponding to error state along with
%visual display before moving onto the specified state.
%   WL_ERROR_STATE(OBJ, MSG, STATE) takes in a pointer to the experiment
%   WLect OBJ displays the error before moving onto state STATE.
    if strcmpi(msg,'Too Fast')
        wl_play_sound(WL, WL.cfg.toofastbeep,WL.cfg.vol);
    elseif ~isempty(msg)
        wl_play_sound(WL, WL.cfg.tooslowbeep,WL.cfg.vol);
    end

    WL.State.ErrorResume = state;
    wl_state_next(WL, WL.State.ERROR);
    WL.GW.error_msg = msg;
end
