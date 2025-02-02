defmodule Wavelets.NDWT do
  @moduledoc """
  Implementation of n-dimensional Discrete Wavelet Transform using Mallat's orthonormal convention
  """

  alias Wavelets.{DWT, Filter, IDWT}

  @doc """
  Performs n-dimensional forward discrete wavelet transform
  """
  def forward(signal, filter, dims, precision \\ :double) when dims > 0 do
    case dims do
      1 when is_list(signal) ->
        DWT.forward_1d(signal, filter, precision)

      2 when is_list(hd(signal)) ->
        DWT.forward_2d(signal, filter, precision)

      _ ->
        signal
        |> Enum.map(&forward(&1, filter, dims - 1, precision))
        |> Enum.unzip()
        |> apply_dimension_scaling(dims)
    end
  end

  @doc """
  Performs n-dimensional inverse discrete wavelet transform
  """
  def inverse(approximation, details, filter, dims, precision \\ :double)
      when dims > 0 do
    # Remove dimension scaling before inverse transform
    {scaled_approx, scaled_details} =
      remove_dimension_scaling(approximation, details, dims)

    case dims do
      1 ->
        IDWT.inverse_1d(scaled_approx, scaled_details, filter, precision)

      2 ->
        IDWT.inverse_2d(scaled_approx, scaled_details, filter, precision)

      _ ->
        Enum.zip(scaled_approx, scaled_details)
        |> Enum.map(fn {a, d} -> inverse(a, d, filter, dims - 1, precision) end)
    end
  end

  # Apply additional scaling for higher dimensions
  defp apply_dimension_scaling({approx, details}, dims) when dims > 2 do
    scale = :math.pow(:math.sqrt(2), dims - 2)

    {
      scale_nd(approx, scale, dims),
      scale_nd(details, scale, dims)
    }
  end

  defp apply_dimension_scaling(result, _), do: result

  # Remove dimension scaling for inverse transform
  defp remove_dimension_scaling(approx, details, dims) when dims > 2 do
    scale = 1 / :math.pow(:math.sqrt(2), dims - 2)

    {
      scale_nd(approx, scale, dims),
      scale_nd(details, scale, dims)
    }
  end

  defp remove_dimension_scaling(approx, details, _), do: {approx, details}

  defp scale_nd(data, scale, dims) do
    case data do
      list when is_list(list) ->
        if dims <= 1 do
          Enum.map(list, &(&1 * scale))
        else
          Enum.map(list, &scale_nd(&1, scale, dims - 1))
        end

      value when is_number(value) ->
        value * scale
    end
  end
end
