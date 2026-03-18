function [tracks, track_props] = filter_duration_short(tracks, track_props, min_duration)

short_tracks = [track_props(:).duration] > min_duration;
tracks(short_tracks) = [];
track_props(short_tracks) = [];
