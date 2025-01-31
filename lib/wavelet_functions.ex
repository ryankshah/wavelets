defmodule Wavelets.WaveletFunctions do
  @moduledoc """
  Module for computing approximations of wavelet and scaling functions
  """

  alias Wavelets.Filter

  @doc """
  Computes points of the scaling function phi(x) for a given wavelet filter
  """
  @spec scaling_function_points(Filter.t(), Range.t()) :: list({number, number})
  def scaling_function_points(filter, x_range) do
    # Initialize with delta function
    init_points =
      Enum.map(x_range, fn x ->
        if x == 0, do: {x, 1.0}, else: {x, 0.0}
      end)

    # Iterate cascade algorithm
    # number of iterations for convergence
    iterations = 7

    refined_points =
      1..iterations
      |> Enum.reduce(init_points, fn _i, points ->
        cascade_algorithm_step(points, filter.reconstruction_low_pass)
      end)

    # Normalize
    sum =
      refined_points
      |> Enum.map(fn {_, y} -> y end)
      |> Enum.sum()

    Enum.map(refined_points, fn {x, y} -> {x, y / sum} end)
  end

  @doc """
  Computes points of the wavelet function psi(x) for a given wavelet filter
  """
  @spec wavelet_function_points(Filter.t(), Range.t()) :: list({number, number})
  def wavelet_function_points(filter, x_range) do
    # First get scaling function points
    phi_points = scaling_function_points(filter, x_range)

    # Compute wavelet function using two-scale relation
    Enum.map(x_range, fn x ->
      y =
        phi_points
        |> Enum.with_index()
        |> Enum.map(fn {{_, phi_y}, k} ->
          coeff = Enum.at(filter.reconstruction_high_pass, k, 0)
          phi_y * coeff
        end)
        |> Enum.sum()

      {x, y}
    end)
  end

  defp cascade_algorithm_step(points, coeffs) do
    dx =
      points
      |> Enum.chunk_every(2, 1)
      |> hd()
      |> then(fn [{x1, _}, {x2, _}] -> x2 - x1 end)

    Enum.map(points, fn {x, _} ->
      y =
        coeffs
        |> Enum.with_index()
        |> Enum.map(fn {c, k} ->
          # Find value of phi(2x - k) by interpolation
          x_shifted = 2 * x - k

          interpolate_point(points, x_shifted)
          |> Kernel.*(c)
        end)
        |> Enum.sum()

      {x, y}
    end)
  end

  defp interpolate_point(points, x) do
    # Linear interpolation between nearest points
    {x1, y1} = Enum.min_by(points, fn {px, _} -> abs(px - x) end)
    {x2, y2} = Enum.min_by(points -- [{x1, y1}], fn {px, _} -> abs(px - x) end)

    if x1 == x2 do
      y1
    else
      y1 + (y2 - y1) * (x - x1) / (x2 - x1)
    end
  end
end
